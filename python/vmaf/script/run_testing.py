#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')

import os
import sys
import re

import numpy as np
from vmaf.config import DisplayConfig

from vmaf.core.result_store import FileSystemResultStore
from vmaf.tools.misc import import_python_file, get_cmd_option, cmd_option_exists
from vmaf.core.quality_runner import QualityRunner, VmafQualityRunner, BootstrapVmafQualityRunner, DlbDeitpVmafQualityRunner
from vmaf.core.matlab_quality_runner import STMADQualityRunner, SpEEDMatlabQualityRunner, StrredQualityRunner, StrredOptQualityRunner
from vmaf.routine import run_test_on_dataset, print_matplotlib_warning
from vmaf.tools.stats import ListStats

__copyright__ = "Copyright 2016-2020, Netflix, Inc."
__license__ = "BSD+Patent"

POOL_METHODS = ['mean', 'harmonic_mean', 'min', 'median', 'perc5', 'perc10', 'perc20']

SUBJECTIVE_MODELS = ['DMOS (default)', 'DMOS_MLE', 'MLE', 'MOS', 'SR_DMOS', 'SR_MOS', 'ZS_SR_DMOS', 'ZS_SR_MOS']


def print_usage():
    quality_runner_types = ['VMAF', 'PSNR', 'SSIM', 'MS_SSIM']
    print("usage: " + os.path.basename(sys.argv[0]) + \
          " quality_type test_dataset_filepath [--vmaf-model VMAF_model_path] " \
          "[--vmaf-phone-model] [--subj-model subjective_model] [--cache-result] " \
          "[--parallelize] [--print-result] [--save-plot plot_dir] [--plot-wh plot_wh]\n")
    print("quality_type:\n\t" + "\n\t".join(quality_runner_types) +"\n")
    print("subjective_model:\n\t" + "\n\t".join(SUBJECTIVE_MODELS) + "\n")
    print("plot_wh: plot width and height in inches, example: 5x5 (default)")


def main():
    if len(sys.argv) < 3:
        print_usage()
        return 2

    try:
        quality_type = sys.argv[1]
        test_dataset_filepath = sys.argv[2]
    except ValueError:
        print_usage()
        return 2

    vmaf_model_path = get_cmd_option(sys.argv, 3, len(sys.argv), '--vmaf-model')
    cache_result = cmd_option_exists(sys.argv, 3, len(sys.argv), '--cache-result')
    parallelize = cmd_option_exists(sys.argv, 3, len(sys.argv), '--parallelize')
    print_result = cmd_option_exists(sys.argv, 3, len(sys.argv), '--print-result')
    suppress_plot = cmd_option_exists(sys.argv, 3, len(sys.argv), '--suppress-plot')
    vmaf_phone_model = cmd_option_exists(sys.argv, 3, len(sys.argv), '--vmaf-phone-model')

    

    pool_method = get_cmd_option(sys.argv, 3, len(sys.argv), '--pool')
    if not (pool_method is None
            or pool_method in POOL_METHODS):
        print('--pool can only have option among {}'.format(', '.join(POOL_METHODS)))
        return 2

    subj_model = get_cmd_option(sys.argv, 3, len(sys.argv), '--subj-model')

    try:
        if subj_model is not None:
            from sureal.subjective_model import SubjectiveModel
            subj_model_class = SubjectiveModel.find_subclass(subj_model)
        else:
            subj_model_class = None
    except Exception as e:
        print("Error: " + str(e))
        return 1

    save_plot_dir = get_cmd_option(sys.argv, 3, len(sys.argv), '--save-plot')

    plot_wh = get_cmd_option(sys.argv, 3, len(sys.argv), '--plot-wh')
    if plot_wh is not None:
        try:
            mo = re.match(r"([0-9]+)x([0-9]+)", plot_wh)
            assert mo is not None
            w = mo.group(1)
            h = mo.group(2)
            w = int(w)
            h = int(h)
            plot_wh = (w, h)
        except Exception as e:
            print("Error: plot_wh must be in the format of WxH, example: 5x5")
            return 1

    try:
        runner_class = QualityRunner.find_subclass(quality_type)
    except Exception as e:
        print("Error: " + str(e))
        return 1

    if vmaf_model_path is not None and runner_class != VmafQualityRunner and \
                        runner_class != BootstrapVmafQualityRunner and \
                        runner_class != DlbDeitpVmafQualityRunner:
        print("Input error: only quality_type of VMAF or DEITP_VMAF accepts --vmaf-model.")
        print_usage()
        return 2

    if vmaf_phone_model and runner_class != VmafQualityRunner and \
                    runner_class != BootstrapVmafQualityRunner:
        print("Input error: only quality_type of VMAF accepts --vmaf-phone-model.")
        print_usage()
        return 2

    try:
        test_dataset = import_python_file(test_dataset_filepath)
    except Exception as e:
        print("Error: " + str(e))
        return 1

    if cache_result:
        result_store = FileSystemResultStore()
    else:
        result_store = None

    # pooling
    if pool_method == 'harmonic_mean':
        aggregate_method = ListStats.harmonic_mean
    elif pool_method == 'min':
        aggregate_method = np.min
    elif pool_method == 'median':
        aggregate_method = np.median
    elif pool_method == 'perc5':
        aggregate_method = ListStats.perc5
    elif pool_method == 'perc10':
        aggregate_method = ListStats.perc10
    elif pool_method == 'perc20':
        aggregate_method = ListStats.perc20
    else: # None or 'mean'
        aggregate_method = np.mean

    if vmaf_phone_model:
        enable_transform_score = True
    else:
        enable_transform_score = None

    try:
        if suppress_plot:
            raise AssertionError

        from vmaf import plt
        if plot_wh is None:
            plot_wh = (5, 5)
        fig, ax = plt.subplots(figsize=plot_wh, nrows=1, ncols=1)

        assets, results = run_test_on_dataset(test_dataset, runner_class, ax,
                                          result_store, vmaf_model_path,
                                          parallelize=parallelize,
                                          aggregate_method=aggregate_method,
                                          subj_model_class=subj_model_class,
                                          enable_transform_score=enable_transform_score
                                          )

        bbox = {'facecolor':'white', 'alpha':0.5, 'pad':20}
        ax.annotate('Testing Set', xy=(0.1, 0.85), xycoords='axes fraction', bbox=bbox)

        # ax.set_xlim([-10, 110])
        # ax.set_ylim([-10, 110])

        plt.tight_layout()

        if save_plot_dir is None:
            DisplayConfig.show()
        else:
           # DisplayConfig.show(write_to_dir=save_plot_dir)
            DisplayConfig.show(write_to_dir=save_plot_dir,
                               data_set_name=test_dataset.dataset_name,
                               model_name=os.path.splitext(vmaf_model_path.split("/")[-1])[0])

    except ImportError:
        print_matplotlib_warning()
        assets, results = run_test_on_dataset(test_dataset, runner_class, None,
                                          result_store, vmaf_model_path,
                                          parallelize=parallelize,
                                          aggregate_method=aggregate_method,
                                          subj_model_class=subj_model_class,
                                          enable_transform_score=enable_transform_score
                                          )
    except AssertionError:
        assets, results = run_test_on_dataset(test_dataset, runner_class, None,
                                          result_store, vmaf_model_path,
                                          parallelize=parallelize,
                                          aggregate_method=aggregate_method,
                                          subj_model_class=subj_model_class,
                                          enable_transform_score=enable_transform_score
                                          )

    save_result_file = get_cmd_option(sys.argv, 3, len(sys.argv), '--save-result-file')
    
    if print_result:
        import pandas as pd
        if save_result_file is not None:
            writer = pd.ExcelWriter(save_result_file, mode='ab')
        for result in results:
            if save_result_file is None:
                pd.set_option("display.max_rows", None,
                              "display.max_columns", None)
                #print(result.to_xml())
                #print(str(result.to_dataframe()))
                #print(str(result.to_dataframe().to_dict()))
                #print('')
                predicted_score = result[runner_class.get_score_key()]
                groundtruth_score = result.asset.groundtruth
                print('asset ID: {} predicted: {} groundtruth: {}'.format(result.asset.asset_id,int(predicted_score),int(groundtruth_score)))
            else:
                list_scores_key = result.get_ordered_list_scores_key()
                list_scores = list(map(lambda key: result.result_dict[key],list_scores_key))
                df = pd.DataFrame()
                for list_score, scores in zip(list_scores_key,list_scores):
                    df[list_score] = scores
                    
                predicted_score = result[runner_class.get_score_key()]
                groundtruth_score = result.asset.groundtruth
                asset_sheet_name = "asset_" + str(result.asset.asset_id) + \
                                    "_ps_" + str(int(predicted_score)) + \
                                    "_gt_" + str(int(groundtruth_score))
                df.to_excel(writer, sheet_name=asset_sheet_name)
                print('asset ID: {} predicted: {} groundtruth: {}'.format(result.asset.asset_id,int(predicted_score),int(groundtruth_score)))
                writer.save()
               # print(result.asset.asset_id)
        if save_result_file is not None:
            if writer:
                writer.close()
      


    return 0


if __name__ == '__main__':
    ret = main()
    exit(ret)
