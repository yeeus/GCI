#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_backend
from matplotlib.ticker import FuncFormatter
import argparse
import os
import gzip
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict
import bisect
from typing import Dict, List, Tuple, Optional
import pandas as pd
from array import array

class BedRegionParser:
    """解析BED格式的区间文件"""

    def __init__(self, bed_file: str):
        self.bed_file = bed_file

    def parse_bed_regions(self) -> Dict[str, List[Tuple[int, int]]]:
        """解析BED文件，返回每个染色体的区间列表"""
        regions = defaultdict(list)

        with open(self.bed_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        chrom = parts[0]
                        start = int(parts[1])
                        end = int(parts[2])
                        regions[chrom].append((start, end))

        # 对每个染色体的区间进行排序
        for chrom in regions:
            regions[chrom].sort()

        return dict(regions)

class DepthParser:
    """解析类似fasta格式的深度文件"""

    def __init__(self, depth_file: str):
        self.depth_file = depth_file
        self.sequences = {}  # 存储 seq_id -> array.array('H')
        self.mean_depths = {}  # 存储每个序列的平均深度
        # self.index 不再需要，因为我们将直接操作 array.array 或其numpy转换

    def parse_depth_file_filtered(self, target_sequences: set) -> Dict[str, array]: # 返回值类型改为 array
        """解析深度文件，只处理target_sequences中指定的序列，使用array.array存储深度"""
        # sequences 字典现在将存储 array.array 对象
        current_seq = None
        # 使用 'H' 表示无符号短整型 (0-65535)。如果深度可能更大或为负，选择其他类型码，如 'L', 'i', 'l'
        # depths_list = [] # 不再使用 python list
        current_depth_array = None
        include_current = False

        opener = gzip.open if self.depth_file.endswith('.gz') else open

        print(f"开始解析深度文件: {self.depth_file}")
        print(f"目标序列数量: {len(target_sequences)}")

        processed_count = 0
        skipped_count = 0

        with opener(self.depth_file, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq and current_depth_array is not None and include_current:
                        self.sequences[current_seq] = current_depth_array
                        # 计算平均深度时需要转换为numpy array或逐个累加
                        # 为了效率，如果数据量大，在需要时再计算，或者在添加完所有数据后计算
                        # 这里我们先存起来，get_mean_depth 时再计算
                        processed_count += 1
                        if processed_count > 0 and processed_count % 5 == 0:
                            print(f"已收集 {processed_count} 个目标序列的深度数据...")
                    elif current_seq and not include_current:
                        skipped_count += 1

                    current_seq = line[1:]
                    include_current = current_seq in target_sequences
                    if include_current:
                        current_depth_array = array('H') # 为新序列创建新的 array
                    else:
                        current_depth_array = None # 如果不包含，则不创建 array

                elif include_current and current_depth_array is not None:
                    try:
                        depth = int(line)
                        # 对于 'H'，确保值在 0-65535 范围内，否则会 OverflowError
                        # 可以添加检查和剪切，或者选择更合适的数据类型
                        if 0 <= depth <= 65535:
                            current_depth_array.append(depth)
                        elif depth > 65535:
                             current_depth_array.append(65535) # Cap at max value for 'H'
                        # else: ignore negative or very small values if not applicable
                    except ValueError:
                        continue

            if current_seq and current_depth_array is not None and include_current:
                self.sequences[current_seq] = current_depth_array
                processed_count += 1
            elif current_seq and not include_current:
                skipped_count += 1

        print(f"解析完成: 收集了 {processed_count} 个序列的深度数据，跳过了 {skipped_count} 个序列")
        # 计算并存储平均深度
        for seq_id, arr in self.sequences.items():
            if len(arr) > 0:
                # np.mean 在 array.array 上直接工作可能效率不高或不支持，转为np array
                self.mean_depths[seq_id] = np.mean(np.array(arr, dtype=np.float32)) # 使用float32计算均值
            else:
                self.mean_depths[seq_id] = 0.0
        return self.sequences # 返回包含 array.array 对象的字典

    # parse_depth_file 方法也应做类似修改，或标记为弃用/内部使用
    # 为了简洁，这里省略 parse_depth_file 的修改，假设主要使用 filtered 版本

    def get_region_depths(self, seq_id: str, start: int, end: int) -> np.ndarray:
        """快速获取特定区间的深度数据，返回NumPy数组"""
        if seq_id not in self.sequences:
            return np.array([])

        depth_array_obj = self.sequences[seq_id]
        # 切片 array.array 对象，然后转换为 NumPy 数组
        # 确保 start 和 end 在有效范围内
        actual_start = max(0, start)
        actual_end = min(len(depth_array_obj), end + 1)
        if actual_start >= actual_end: # 如果区间无效或为空
            return np.array([])

        return np.array(depth_array_obj[actual_start:actual_end], dtype=np.int32)

    def get_mean_depth(self, seq_id: str) -> float:
        """获取序列的平均深度"""
        return self.mean_depths.get(seq_id, 0.0)

class FAIParser:
    """解析fai索引文件"""

    def __init__(self, fai_file: str):
        self.fai_file = fai_file

    def parse_fai(self) -> Dict[str, int]:
        """解析fai文件，返回序列ID到长度的映射"""
        seq_lengths = {}

        with open(self.fai_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    seq_id = parts[0]
                    length = int(parts[1])
                    seq_lengths[seq_id] = length

        return seq_lengths

class SlidingWindowProcessor:
    """滑动窗口处理和合并"""

    def __init__(self, window_size: int = 1000):
        self.window_size = window_size

    def calculate_sliding_window(self, depths_input) -> Tuple[np.ndarray, np.ndarray]: # depths_input可以是array.array或np.ndarray
        """计算滑动窗口平均深度"""

        # 如果输入是 array.array，先转换为 NumPy array
        if isinstance(depths_input, array):
            depths = np.array(depths_input, dtype=np.int32) # 或者根据需要选择更合适dtype
        elif isinstance(depths_input, np.ndarray):
            depths = depths_input
        else:
            raise TypeError("Input depths must be an array.array or numpy.ndarray")

        if len(depths) == 0:
            return np.array([]), np.array([])

        # 使用numpy的卷积进行快速滑动窗口计算
        kernel = np.ones(self.window_size) / self.window_size

        # 处理边界情况
        if len(depths) < self.window_size:
            # 如果深度数据长度小于窗口大小，则整个区域的平均深度作为单个点
            avg_depths = np.array([np.mean(depths)])
            # 位置设为区域中心
            positions = np.array([len(depths) // 2]) 
        else:
            # 使用valid模式避免边界效应
            avg_depths = np.convolve(depths, kernel, mode='valid')
            # 计算窗口中心位置
            # 第一个窗口的中心是 window_size // 2
            # 最后一个窗口的中心是 len(depths) - (window_size // 2) -1 (如果window_size是奇数)
            # 或者 len(depths) - (window_size / 2) (如果window_size是偶数)
            # np.convolve(..., mode='valid') 的输出长度是 M - N + 1，其中M是depths长度，N是window_size
            # 我们需要为这 M - N + 1 个点找到中心位置
            # 第一个有效窗口覆盖 depths[0:window_size]，其中心应为 (0 + window_size - 1) // 2 = (window_size - 1) // 2
            # 最后一个有效窗口覆盖 depths[len(depths)-window_size : len(depths)]，其中心应为 (len(depths)-window_size + len(depths)-1)//2

            # 修正 positions 的计算
            # 第一个窗口的起始索引是0，结束索引是 window_size - 1。中心点是 (window_size - 1) // 2。
            # 最后一个窗口的起始索引是 len(depths) - window_size，结束索引是 len(depths) - 1。
            # 中心点是 (len(depths) - window_size + len(depths) - 1) // 2。
            # 生成的 avg_depths 长度为 len(depths) - self.window_size + 1
            # positions 应该与 avg_depths 的长度相同
            start_pos = (self.window_size - 1) // 2
            # 最后一个有效窗口的中心位置
            # 最后一个窗口的起始点是 len(depths) - self.window_size
            # 它的中心点是 (len(depths) - self.window_size) + (self.window_size - 1) // 2
            end_pos = (len(depths) - self.window_size) + (self.window_size - 1) // 2

            # 如果 avg_depths 为空 (例如 window_size > len(depths) 但不满足 len(depths) < self.window_size 的初始条件)
            if len(avg_depths) == 0:
                # 这种情况理论上不应该发生，因为上面已经有 len(depths) < self.window_size 的判断
                # 但为了安全起见
                if len(depths) > 0:
                    avg_depths = np.array([np.mean(depths)])
                    positions = np.array([len(depths) // 2])
                else:
                    return np.array([]), np.array([]) # 如果depths为空，直接返回空
            else:
                 # positions 应该是每个有效窗口的中心点索引
                 # 第一个有效窗口是 depths[0:window_size], 中心是 (window_size-1)//2
                 # 第二个有效窗口是 depths[1:window_size+1], 中心是 1 + (window_size-1)//2
                 # ...
                 # 最后一个有效窗口是 depths[L-W : L], 中心是 (L-W) + (window_size-1)//2
                 # L = len(depths), W = window_size
                num_valid_windows = len(depths) - self.window_size + 1
                positions = np.arange(num_valid_windows) + (self.window_size -1) // 2

        return positions, avg_depths

    def merge_consecutive_windows(self, positions: np.ndarray,
                                depths: np.ndarray) -> Tuple[List[Tuple[int, int]], List[float]]:
        """合并连续相同深度的窗口"""
        if len(depths) == 0:
            return [], []

        merged_regions = []
        merged_depths = []

        current_start = positions[0]
        current_depth = depths[0]
        current_end = positions[0]

        for i in range(1, len(depths)):
            if abs(depths[i] - current_depth) < 0.1:  # 允许小的浮点误差
                current_end = positions[i]
            else:
                # 保存当前区间
                merged_regions.append((current_start, current_end))
                merged_depths.append(current_depth)

                # 开始新区间
                current_start = positions[i]
                current_depth = depths[i]
                current_end = positions[i]

        # 保存最后一个区间
        merged_regions.append((current_start, current_end))
        merged_depths.append(current_depth)

        return merged_regions, merged_depths

class DepthRegionAnalyzer:
    """Analyze depth regions to identify problematic areas"""

    def __init__(self, min_safe_depth: int = 5):
        self.min_safe_depth = min_safe_depth

    def analyze_depth_regions(self, hifi_depths: np.ndarray, ont_depths: np.ndarray):
        """Analyze depth data to identify zero-depth, low-depth, and normal regions"""
        # Combine HiFi and ONT depths for analysis
        if len(ont_depths) > 0 and len(hifi_depths) > 0:
            combined_depths = hifi_depths + ont_depths
        elif len(ont_depths) > 0:
            combined_depths = ont_depths
        else:
            combined_depths = hifi_depths

        # Find different types of regions
        zero_depth_regions = self._find_continuous_regions(combined_depths == 0)
        low_depth_regions = self._find_continuous_regions(
            (combined_depths > 0) & (combined_depths < self.min_safe_depth)
        )
        normal_depth_regions = self._find_continuous_regions(combined_depths >= self.min_safe_depth)

        return {
            'zero_depth': zero_depth_regions,
            'low_depth': low_depth_regions,
            'normal_depth': normal_depth_regions
        }

    def _find_continuous_regions(self, mask: np.ndarray):
        """Find continuous regions where mask is True"""
        if len(mask) == 0:
            return []

        regions = []
        start = None

        for i, value in enumerate(mask):
            if value and start is None:
                start = i
            elif not value and start is not None:
                regions.append((start, i - 1))
                start = None

        # Handle case where region extends to the end
        if start is not None:
            regions.append((start, len(mask) - 1))

        return regions

class DataProcessor:
    """统一的数据处理器，处理单一类型的深度数据"""

    def __init__(self, data_type: str, color: str, window_size: int = 1000):
        self.data_type = data_type  # 'hifi' or 'ont'
        self.color = color
        self.window_size = window_size

    def calculate_windowed_stats(self, depths: np.ndarray):
        """计算窗口统计信息"""
        if len(depths) == 0:
            return np.array([]), np.array([]), np.array([])

        # 找到非零区域段
        zero_mask = (depths == 0)
        segments = self._find_non_zero_segments(depths, zero_mask)

        all_means = []
        all_starts = []
        all_ends = []

        for segment_start, segment_end in segments:
            segment_depths = depths[segment_start:segment_end+1]
            segment_length = len(segment_depths)

            if segment_length == 0:
                continue

            num_windows = max(1, segment_length // self.window_size)
            if segment_length % self.window_size != 0:
                num_windows += 1

            for i in range(num_windows):
                window_start_in_segment = i * self.window_size
                window_end_in_segment = min((i + 1) * self.window_size, segment_length)

                abs_start = segment_start + window_start_in_segment
                abs_end = segment_start + window_end_in_segment - 1

                window_data = segment_depths[window_start_in_segment:window_end_in_segment]
                if len(window_data) > 0:
                    all_means.append(np.mean(window_data))
                    all_starts.append(abs_start)
                    all_ends.append(abs_end)

        return np.array(all_means), np.array(all_starts), np.array(all_ends)

    def _find_non_zero_segments(self, depths: np.ndarray, zero_mask: np.ndarray):
        """找到连续的非零段"""
        segments = []
        start = None

        for i, is_zero in enumerate(zero_mask):
            if not is_zero and start is None:
                start = i
            elif is_zero and start is not None:
                segments.append((start, i - 1))
                start = None

        if start is not None:
            segments.append((start, len(depths) - 1))

        return segments

    def analyze_depth_regions(self, depths: np.ndarray, min_safe_depth: int = 5):
        """分析深度区域（独立计算），包括深度为0的区域"""
        if len(depths) == 0:
            return {'zero': [], 'low': [], 'medium': []}

        # 定义深度阈值
        low_threshold = min_safe_depth

        # 创建掩码
        zero_mask = depths == 0
        low_mask = (depths > 0) & (depths < low_threshold)

        return {
            'zero': self._mask_to_regions(zero_mask),
            'low': self._mask_to_regions(low_mask),
        }

    def _mask_to_regions(self, mask: np.ndarray):
        """将布尔掩码转换为区域列表"""
        regions = []
        start = None

        for i, value in enumerate(mask):
            if value and start is None:
                start = i
            elif not value and start is not None:
                regions.append((start, i - 1))
                start = None

        if start is not None:
            regions.append((start, len(mask) - 1))

        return regions

class DepthPlotter:
    def __init__(self, hifi_color: str = '#2ca25f', ont_color: str = '#3C5488',
                 figure_size: tuple = (15, 4), output_format: str = 'png', dpi: int = 300,
                 max_depth_ratio: float = 4.0):
        self.hifi_color = hifi_color
        self.ont_color = ont_color
        self.figure_size = figure_size
        self.output_format = output_format
        self.dpi = dpi
        self.output_dir = '.'
        self.max_depth_ratio = max_depth_ratio

    def _setup_plot_properties(self, ax, seq_id: str, seq_length: int, plot_mode: str, avg_depth: float = None):
        """set basic plot properties"""
        # set title
        ax.set_title(f'Depth Coverage for {seq_id}', fontsize=14, fontweight='bold')

        # set X axis
        formatter = FuncFormatter(self._format_position_axis)
        ax.xaxis.set_major_formatter(formatter)
        position_label = self._get_position_unit_label(seq_length)
        ax.set_xlabel(position_label, fontsize=12)
        ax.set_xlim(0, seq_length)

        # set Y axis
        ax.set_ylabel('Depth', fontsize=12)

        # set Y range according to data type
        if plot_mode == 'both':
            # double Y, positive and negative
            if avg_depth is not None and hasattr(self, 'max_depth_ratio'):
                max_y = avg_depth * self.max_depth_ratio
                ax.set_ylim(-max_y, max_y)
            else:
                ax.set_ylim(-100, 100)
            # y=0
            ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.7)

            # 使用自定义格式化器来消除警告：
            def abs_formatter(x, pos):
                return str(abs(int(x)))
            ax.yaxis.set_major_formatter(FuncFormatter(abs_formatter))
        else:
            # y is positive
            if avg_depth is not None and hasattr(self, 'max_depth_ratio'):
                max_y = avg_depth * self.max_depth_ratio
                ax.set_ylim(0, max_y)
            else:
                ax.set_ylim(0, 100)

        # add grid
        ax.grid(True, alpha=0.2)


    def plot_single_sequence(self, seq_id: str, hifi_depths: List[int], ont_depths: List[int],
                           window_size: int = 1000, regions: list = None, output_dir: str = None,
                           min_safe_depth: int = 5) -> str:
        """plot depth plot for single sequence, support 3 models"""

        # convert to numpy
        hifi_depths_array = np.array(hifi_depths) if hifi_depths is not None else np.array([])
        ont_depths_array = np.array(ont_depths) if ont_depths is not None else np.array([])

        # check data valid
        has_hifi = len(hifi_depths_array) > 0
        has_ont = len(ont_depths_array) > 0

        if not has_hifi and not has_ont:
            print(f"Error: No depth data for sequence {seq_id}")
            return None

        # check length coordance when both exists
        if has_hifi and has_ont:
            if len(hifi_depths_array) != len(ont_depths_array):
                raise ValueError(f"Error: HiFi and ONT data length mismatch for sequence {seq_id}. "
                               f"HiFi length: {len(hifi_depths_array)}, ONT length: {len(ont_depths_array)}. "
                               f"Both datasets must have the same length.")

        # sequence length
        seq_length = max(len(hifi_depths_array), len(ont_depths_array))

        # creat a data processor
        processors = {}
        if has_hifi:
            processors['hifi'] = DataProcessor('hifi', self.hifi_color, window_size)
        if has_ont:
            processors['ont'] = DataProcessor('ont', self.ont_color, window_size)

        # processing data
        processed_data = {}
        for data_type, processor in processors.items():
            if data_type == 'hifi':
                depths = hifi_depths_array
            else:
                depths = ont_depths_array

            # calculate depth by windows
            means, starts, ends = processor.calculate_windowed_stats(depths)

            # analyze depth categories
            depth_regions = processor.analyze_depth_regions(depths, min_safe_depth)

            processed_data[data_type] = {
                'means': means,
                'starts': starts,
                'ends': ends,
                'regions': depth_regions,
                'processor': processor
            }

        # set plot mode
        if has_hifi and has_ont:
            plot_mode = 'both'  # HiFi is upper, ONT is lower
        elif has_hifi:
            plot_mode = 'hifi_only'  # only HiFi, Y is positive
        else:
            plot_mode = 'ont_only'  # only ONT, Y is positive

        # create a figure
        fig, ax = plt.subplots(figsize=self.figure_size, dpi=self.dpi)

        # ploting based on mode
        if plot_mode == 'both':
            self._plot_both_data(ax, processed_data, seq_length)
        elif plot_mode == 'hifi_only':
            self._plot_single_data(ax, processed_data['hifi'], seq_length)
        else:  # ont_only
            self._plot_single_data(ax, processed_data['ont'], seq_length)

        # calculating the average depth for considering Y range
        all_depths = []
        if has_hifi:
            all_depths.extend(hifi_depths_array[hifi_depths_array > 0])  # exculde reigons with 0 depth
        if has_ont:
            all_depths.extend(ont_depths_array[ont_depths_array > 0])  # exculde reigons with 0 depth

        avg_depth = np.mean(all_depths) if len(all_depths) > 0 else 1.0

        # pass average depth information to plot properties
        self._setup_plot_properties(ax, seq_id, seq_length, plot_mode, avg_depth)

        # add legend
        self._create_legend(ax, processed_data, plot_mode)

        # save plot
        output_path = self._save_plot(fig, seq_id, regions, output_dir)

        # close plot
        plt.close(fig)

        return output_path

    def _plot_both_data(self, ax, processed_data, seq_length):
        """plot hifi and ont data（HiFi upper，ONT lower）"""
        # for hifi, positive
        if 'hifi' in processed_data:
            data = processed_data['hifi']
            self._plot_depth_regions(ax, data['regions'], positive=True)
            self._plot_windowed_data(ax, data, positive=True)

        # for ont, negative
        if 'ont' in processed_data:
            data = processed_data['ont']
            self._plot_depth_regions(ax, data['regions'], positive=False)
            self._plot_windowed_data(ax, data, positive=False)

    def _plot_single_data(self, ax, data, seq_length):
        """plot one dateset"""
        self._plot_depth_regions(ax, data['regions'], positive=True)
        self._plot_windowed_data(ax, data, positive=True)

    def _plot_windowed_data(self, ax, data, positive=True):
        """plot windows data"""
        means = data['means']
        starts = data['starts']
        ends = data['ends']
        color = data['processor'].color

        if len(means) == 0:
            return

        # calculate the middle position of any window
        centers = (starts + ends) / 2

        # plot upper or lower
        plot_means = means if positive else -means

        # plot bar 
        widths = ends - starts + 1
        ax.bar(centers, plot_means, width=widths, color=color, alpha=0.8, edgecolor='none')

        # add average line
        if len(means) > 0:
            avg_depth = np.mean(means)
            avg_line = avg_depth if positive else -avg_depth
            ax.axhline(y=avg_line, color=color, linestyle='--', alpha=0.8, linewidth=1)

    def _plot_depth_regions(self, ax, regions, positive=True):
        """plot background for different kinds of depth regions"""
        region_colors = {
            'zero': '#FAD7DD',      # depth = 0
            'low': '#B7DBEA',       # low depth
        }

        y_min, y_max = ax.get_ylim()

        for region_type, region_list in regions.items():
            if region_type not in region_colors:
                continue

            color_value = region_colors[region_type]

            for start, end in region_list:
                if positive:
                    ax.axvspan(start, end, ymin=0.5, ymax=0.95,
                               color=color_value, alpha=0.8)
                else:
                    ax.axvspan(start, end, ymin=0.05, ymax=0.5,
                               color=color_value, alpha=0.8)

    def _create_legend(self, ax, processed_data, plot_mode):
        """creat legend"""
        legend_elements = []

        for data_type, data in processed_data.items():
            color = data['processor'].color
            label = 'HiFi' if data_type == 'hifi' else 'ONT'
            legend_elements.append(plt.Rectangle((0, 0), 1, 1, facecolor=color, alpha=0.8, label=label))

        # add legend for depth categories, including zero depth
        legend_elements.extend([
            plt.Rectangle((0, 0), 1, 1, facecolor='#FAD7DD', alpha=1.0, label='Zero Depth'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#B7DBEA', alpha=0.8, label='Low Depth')
        ])

        # put legend on upper center
        ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 0.98),
            ncol=len(legend_elements), frameon=True, fancybox=False, shadow=False)

    def _format_position_axis(self, x, pos):
        """normlize x axis"""
        if x >= 1000000:
            return f'{x/1000000:.1f}M'
        elif x >= 1000:
            return f'{x/1000:.1f}k'
        else:
            return f'{int(x)}'

    def _get_position_unit_label(self, max_position):
        """make unit label for x axis"""
        if max_position >= 1000000:
            return 'Position (Mbp)'
        elif max_position >= 1000:
            return 'Position (kbp)'
        else:
            return 'Position (bp)'

    def _save_plot(self, fig, seq_id, regions, output_dir):
        """save figure"""
        if output_dir:
            self.output_dir = output_dir

        if regions:
            region_str = f"_{regions[0][0]}-{regions[0][1]}"
            filename = f"{seq_id}{region_str}.{self.output_format}"
        else:
            filename = f"{seq_id}.{self.output_format}"

        output_path = os.path.join(self.output_dir, filename)
        fig.savefig(output_path, dpi=self.dpi, bbox_inches='tight', facecolor='white')

        return output_path


class SynchronizedDepthReader:
    """Iterator for synchronized reading of two depth files"""

    def __init__(self, hifi_file: str = None, ont_file: str = None, target_sequences: set = None, regions: dict = None):
        self.hifi_file = hifi_file
        self.ont_file = ont_file
        self.target_sequences = target_sequences or set()
        self.regions = regions
        # Add tracking of processed sequences
        self.processed_sequences = set()
        self.total_target_sequences = len(self.target_sequences)

    def read_sequences(self):
        """Synchronously read two files, yield data for each sequence"""
        hifi_opener = gzip.open if self.hifi_file and self.hifi_file.endswith('.gz') else open
        ont_opener = gzip.open if self.ont_file and self.ont_file.endswith('.gz') else open

        hifi_file = hifi_opener(self.hifi_file, 'rt') if self.hifi_file else None
        ont_file = ont_opener(self.ont_file, 'rt') if self.ont_file else None

        try:
            current_seq = None
            hifi_depths = []
            ont_depths = []

            # Read both files simultaneously
            hifi_lines = iter(hifi_file) if hifi_file else iter([])
            ont_lines = iter(ont_file) if ont_file else iter([])

            while True:
                # Check if all target sequences have been processed
                if self.target_sequences and len(self.processed_sequences) >= len(self.target_sequences):
                    print(f"All {len(self.target_sequences)} target sequences have been processed, ending file reading early")
                    break

                try:
                    hifi_line = next(hifi_lines) if hifi_file else None
                    ont_line = next(ont_lines) if ont_file else None

                    if hifi_line is None and ont_line is None:
                        break

                    # Process sequence headers
                    new_seq_from_hifi = None
                    new_seq_from_ont = None

                    # Check for HiFi sequence header
                    if hifi_line and hifi_line.strip().startswith('>'):
                        new_seq_from_hifi = hifi_line.strip()[1:]

                    # Check for ONT sequence header
                    if ont_line and ont_line.strip().startswith('>'):
                        new_seq_from_ont = ont_line.strip()[1:]

                    # If we found a new sequence header
                    if new_seq_from_hifi or new_seq_from_ont:
                        # Process previous sequence first
                        if current_seq and self._should_process_sequence(current_seq):
                            self.processed_sequences.add(current_seq)
                            remaining = len(self.target_sequences) - len(self.processed_sequences)
                            print(f"Processing sequence: {current_seq}, remaining target sequences: {remaining}")
                            yield current_seq, np.array(hifi_depths), np.array(ont_depths)

                            # Check if all target sequences have been processed
                            if self.target_sequences and len(self.processed_sequences) >= len(self.target_sequences):
                                print(f"All target sequences have been processed, stopping reading")
                                break

                        # Set new sequence (prefer HiFi if both are present)
                        current_seq = new_seq_from_hifi or new_seq_from_ont
                        hifi_depths = []
                        ont_depths = []
                    else:
                        # Process depth data
                        if hifi_line and not hifi_line.strip().startswith('>'):
                            try:
                                hifi_depths.append(int(hifi_line.strip()))
                            except ValueError:
                                hifi_depths.append(0)

                        if ont_line and not ont_line.strip().startswith('>'):
                            try:
                                ont_depths.append(int(ont_line.strip()))
                            except ValueError:
                                ont_depths.append(0)

                except StopIteration:
                    break

            # Process the last sequence
            if current_seq and self._should_process_sequence(current_seq) and current_seq not in self.processed_sequences:
                self.processed_sequences.add(current_seq)
                print(f"Processing last sequence: {current_seq}")
                yield current_seq, np.array(hifi_depths), np.array(ont_depths)

        finally:
            if hifi_file:
                hifi_file.close()
            if ont_file:
                ont_file.close()

            print(f"File reading ended, processed {len(self.processed_sequences)} sequences in total")

    def _should_process_sequence(self, seq_id: str) -> bool:
        """Determine whether this sequence should be processed"""
        if self.target_sequences and seq_id not in self.target_sequences:
            return False
        if self.regions and seq_id not in self.regions:
            return False
        return True


def main():
    # Argument parsing
    import argparse
    parser = argparse.ArgumentParser(
        description='Depth data visualization tool - Enhanced version',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Required parameters
    parser.add_argument('-r', '--fai', required=True,
                       help='Reference genome fai index file')

    # Depth file parameters
    parser.add_argument('--hifi',
                       help='HiFi depth file (supports .gz compression)')
    parser.add_argument('--nano',
                       help='ONT depth file (supports .gz compression)')

    # Region parameters
    parser.add_argument('--regions',
                       help='BED format region file')
    parser.add_argument('--region',
                       help='Single region, format: chr:start-end')

    # Output parameters
    parser.add_argument('-o', '--output_dir', default='images',
                       help='Output directory (default: images)')
    parser.add_argument('-f', '--output-format', choices=['png', 'pdf', 'svg'], default='pdf',
                       help='Output format (default: pdf)')

    # Processing parameters
    parser.add_argument('-w', '--window-size', type=int, default=1000,
                       help='Sliding window size (default: 1000)')

    # Depth range parameters
    parser.add_argument('--max-depth-ratio', type=float, default=3.0,
                       help='Maximum depth ratio (relative to average depth, default: 3.0)')

    # New parameters
    parser.add_argument('--min-safe-depth', type=int, default=5,
                       help='Minimum safe depth threshold, regions below this value will be marked with blue background (default: 5)')

    args = parser.parse_args()

    # Check input parameters
    if not args.hifi and not args.nano:
        print("Error: Must provide at least one depth file (--hifi or --nano)")
        return

    # Create output directory
    import os
    os.makedirs(args.output_dir, exist_ok=True)

    # Parse fai file
    print("Parsing fai file...")
    fai_parser = FAIParser(args.fai)
    fai_lengths = fai_parser.parse_fai()
    target_sequences = set(fai_lengths.keys())
    print(f"Found {len(fai_lengths)} reference sequences")

    # Parse region file
    bed_regions = None
    if args.regions:
        print(f"Parsing BED region file: {args.regions}")
        bed_parser = BedRegionParser(args.regions)
        bed_regions = bed_parser.parse_bed_regions()
        total_regions = sum(len(regions) for regions in bed_regions.values())
        print(f"Found {total_regions} regions, involving {len(bed_regions)} sequences")

    # Parse single region parameter
    single_region = None
    if args.region:
        try:
            parts = args.region.split(':')
            seq_id = parts[0]
            start, end = map(int, parts[1].split('-'))
            single_region = {seq_id: [(start, end)]}
            target_sequences = {seq_id}
            print(f"Will plot single specified region: {args.region}")
        except:
            print(f"Error: Invalid region format {args.region}")
            return

    # Prioritize single region, otherwise use BED file regions
    regions_to_use = single_region if single_region else bed_regions

    # Determine target sequences
    if regions_to_use:
        target_sequences = set(regions_to_use.keys())
        print(f"Will process specified regions of {len(target_sequences)} sequences")
    else:
        target_sequences = set(fai_lengths.keys())
        print(f"Will process all {len(target_sequences)} reference sequences")

    # Create synchronized reader
    reader = SynchronizedDepthReader(
        hifi_file=args.hifi,
        ont_file=args.nano,
        target_sequences=target_sequences,
        regions=regions_to_use
    )

    # Create plotter with correct parameters
    plotter = DepthPlotter(output_format=args.output_format)
    plotter.output_dir = args.output_dir

    print("Starting sequential processing...")
    successful = 0
    failed = 0

    # Process each sequence sequentially
    for seq_id, hifi_depths, ont_depths in reader.read_sequences():
        try:
            print(f"Processing sequence: {seq_id}")

            # Determine sequence length from available data
            seq_length = 0
            if len(hifi_depths) > 0:
                seq_length = len(hifi_depths)
            elif len(ont_depths) > 0:
                seq_length = len(ont_depths)
            else:
                print(f"Warning: No depth data for sequence {seq_id}")
                failed += 1
                continue

            # Determine regions to plot
            if regions_to_use and seq_id in regions_to_use:
                sequence_regions = regions_to_use[seq_id]
            else:
                # Plot entire sequence
                sequence_regions = [(0, seq_length - 1)]

            # Generate images for each region
            for region_start, region_end in sequence_regions:
                # Ensure region bounds are valid
                region_start = max(0, region_start)
                region_end = min(seq_length - 1, region_end)

                if region_start > region_end:
                    print(f"Warning: Invalid region [{region_start}, {region_end}] for sequence {seq_id}")
                    continue

                # Extract region data safely
                region_hifi = hifi_depths[region_start:region_end+1] if len(hifi_depths) > 0 else []
                region_ont = ont_depths[region_start:region_end+1] if len(ont_depths) > 0 else []

                result = plotter.plot_single_sequence(
                    seq_id=seq_id,
                    hifi_depths=region_hifi,
                    ont_depths=region_ont,
                    window_size=args.window_size,
                    regions=[(region_start, region_end)],
                    output_dir=args.output_dir
                )

                if result:
                    successful += 1
                    print(f"  Generated: {result}")
                else:
                    failed += 1

        except Exception as e:
            print(f"Error processing sequence {seq_id}: {e}")
            failed += 1

    print(f"\nProcessing completed!")
    print(f"Successful: {successful}, Failed: {failed}")


if __name__ == '__main__':
    main()
