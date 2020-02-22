try:
    from scipy import stats, inf, log
except ImportError:
    class stats:
        class norm:
            def logpdf(foo1, foo2, foo3):
                raise ImportError("No module named scipy")
    inf = float("inf")
    from math import log
from copy import deepcopy
try:
    import numpy as np
except ImportError:
    class np:
        def average(foo1, weights):
            raise ImportError("No module named numpy")
        def arrary(foo1):
            raise ImportError("No module named numpy")
        def vstack(foo1):
            raise ImportError("No module named numpy")
        def where(foo1):
            raise ImportError("No module named numpy")
import random
import sys


def bic(loglike, len_param, len_data):
    return log(len_data) * len_param - 2 * loglike


def aic(loglike, len_param):
    return 2 * len_param - 2 * loglike


def weighted_mean_and_std(values, weights):
    mean = np.average(values, weights=weights)
    std = np.average((values-mean)**2, weights=weights)**0.5
    return mean, std


def weighted_gmm_with_em_aic(data_array, data_weights=None, minimum_cluster=1, maximum_cluster=5, min_sigma_factor=1E-5,
                             cluster_limited=None, log_handler=None, verbose_log=False):
    """
    :param data_array:
    :param data_weights:
    :param minimum_cluster:
    :param maximum_cluster:
    :param min_sigma_factor:
    :param cluster_limited: {dat_id1: {0, 1}, dat_id2: {0}, dat_id3: {0} ...}
    :param log_handler:
    :param verbose_log:
    :return:
    """
    min_sigma = min_sigma_factor * np.average(data_array, weights=data_weights)

    def model_loglike(dat_arr, dat_w, lbs, parameters):
        total_loglike = 0
        for go_to_cl, pr in enumerate(parameters):
            points = dat_arr[lbs == go_to_cl]
            weights = dat_w[lbs == go_to_cl]
            if len(points):
                total_loglike += sum(stats.norm.logpdf(points, pr["mu"], pr["sigma"]) * weights + log(pr["percent"]))
        return total_loglike

    def assign_cluster_labels(dat_arr, dat_w, parameters, limited):
        # assign every data point to its most likely cluster
        if len(parameters) == 1:
            return np.array([0] * int(data_len))
        else:
            # the parameter set of the first cluster
            loglike_res = stats.norm.logpdf(dat_arr, parameters[0]["mu"], parameters[1]["sigma"]) * dat_w + \
                          log(parameters[1]["percent"])
            # the parameter set of the rest cluster
            for pr in parameters[1:]:
                loglike_res = np.vstack(
                    (loglike_res, stats.norm.logpdf(dat_arr, pr["mu"], pr["sigma"]) * dat_w + log(pr["percent"])))
            # assign labels
            new_labels = loglike_res.argmax(axis=0)
            if limited:
                intermediate_labels = []
                for here_dat_id in range(int(data_len)):
                    if here_dat_id in limited:
                        if new_labels[here_dat_id] in limited[here_dat_id]:
                            intermediate_labels.append(new_labels[here_dat_id])
                        else:
                            intermediate_labels.append(sorted(limited[here_dat_id])[0])
                    else:
                        intermediate_labels.append(new_labels[here_dat_id])
                new_labels = np.array(intermediate_labels)
                # new_labels = np.array([
                # sorted(cluster_limited[dat_item])[0]
                # if new_labels[here_dat_id] not in cluster_limited[dat_item] else new_labels[here_dat_id]
                # if dat_item in cluster_limited else
                # new_labels[here_dat_id]
                # for here_dat_id, dat_item in enumerate(data_array)])
                limited_values = set(dat_arr[list(limited)])
            else:
                limited_values = set()
            # re-pick if some cluster are empty
            label_counts = {lb: 0 for lb in range(len(parameters))}
            for ct_lb in new_labels:
                label_counts[ct_lb] += 1
            for empty_lb in label_counts:
                if label_counts[empty_lb] == 0:
                    affordable_lbs = {af_lb: [min, max] for af_lb in label_counts if label_counts[af_lb] > 1}
                    for af_lb in sorted(affordable_lbs):
                        these_points = dat_arr[new_labels == af_lb]
                        if max(these_points) in limited_values:
                            affordable_lbs[af_lb].remove(max)
                        if min(these_points) in limited_values:
                            affordable_lbs[af_lb].remove(min)
                        if not affordable_lbs[af_lb]:
                            del affordable_lbs[af_lb]
                    if affordable_lbs:
                        chose_lb = random.choice(list(affordable_lbs))
                        chose_points = dat_arr[new_labels == chose_lb]
                        data_point = random.choice(affordable_lbs[chose_lb])(chose_points)
                        transfer_index = np.where(dat_arr == data_point)[0]
                        new_labels[transfer_index] = empty_lb
                        label_counts[chose_lb] -= len(transfer_index)
            return new_labels

    def updating_parameter(dat_arr, dat_w, lbs, parameters):

        for go_to_cl, pr in enumerate(parameters):
            these_points = dat_arr[lbs == go_to_cl]
            these_weights = dat_w[lbs == go_to_cl]
            if len(these_points) > 1:
                this_mean, this_std = weighted_mean_and_std(these_points, these_weights)
                pr["mu"] = this_mean
                pr["sigma"] = max(this_std, min_sigma)
                pr["percent"] = sum(these_weights)  # / data_len
            elif len(these_points) == 1:
                pr["sigma"] = max(dat_arr.std() / data_len, min_sigma)
                pr["mu"] = np.average(these_points, weights=these_weights) + pr["sigma"] * (2 * random.random() - 1)
                pr["percent"] = sum(these_weights)  # / data_len
            else:
                # exclude
                pr["mu"] = max(dat_arr) * 1E4
                pr["sigma"] = min_sigma
                pr["percent"] = 1E-10
        return parameters

    data_array = np.array(data_array)
    data_len = float(len(data_array))
    if not len(data_weights):
        data_weights = np.array([1. for foo in range(int(data_len))])
    else:
        assert len(data_weights) == data_len
        average_weights = float(sum(data_weights)) / data_len
        # normalized
        data_weights = np.array([raw_w / average_weights for raw_w in data_weights])

    results = []
    if cluster_limited:
        cls = set()
        for sub_cls in cluster_limited.values():
            cls |= sub_cls
        freedom_dat_item = int(data_len) - len(cluster_limited) + len(cls)
    else:
        freedom_dat_item = int(data_len)
    minimum_cluster = min(freedom_dat_item, minimum_cluster)
    maximum_cluster = min(freedom_dat_item, maximum_cluster)
    for total_cluster_num in range(minimum_cluster, maximum_cluster + 1):
        # initialization
        labels = np.random.choice(total_cluster_num, int(data_len))
        if cluster_limited:
            temp_labels = []
            for dat_id in range(int(data_len)):
                if dat_id in cluster_limited:
                    if labels[dat_id] in cluster_limited[dat_id]:
                        temp_labels.append(labels[dat_id])
                    else:
                        temp_labels.append(sorted(cluster_limited[dat_id])[0])
                else:
                    temp_labels.append(labels[dat_id])
            labels = np.array(temp_labels)
        norm_parameters = updating_parameter(data_array, data_weights, labels,
                                             [{"mu": 0, "sigma": 1, "percent": total_cluster_num/data_len}
                                              for foo in range(total_cluster_num)])
        loglike_shift = inf
        prev_loglike = -inf
        epsilon = 0.01
        count_iterations = 0
        best_loglike = prev_loglike
        best_parameter = norm_parameters
        try:
            while loglike_shift > epsilon:
                count_iterations += 1
                # expectation
                labels = assign_cluster_labels(data_array, data_weights, norm_parameters, cluster_limited)
                # maximization
                updated_parameters = updating_parameter(data_array, data_weights, labels, deepcopy(norm_parameters))
                # loglike shift
                this_loglike = model_loglike(data_array, data_weights, labels, updated_parameters)
                loglike_shift = abs(this_loglike - prev_loglike)
                # update
                prev_loglike = this_loglike
                norm_parameters = updated_parameters
                if this_loglike > best_loglike:
                    best_parameter = updated_parameters
                    best_loglike = this_loglike
            labels = assign_cluster_labels(data_array, data_weights, best_parameter, None)
            results.append({"loglike": best_loglike, "iterates": count_iterations, "cluster_num": total_cluster_num,
                            "parameters": best_parameter, "labels": labels,
                            "aic": aic(prev_loglike, 2 * total_cluster_num),
                            "bic": bic(prev_loglike, 2 * total_cluster_num, data_len)})
        except TypeError as e:
            if log_handler:
                log_handler.error("This error might be caused by outdated version of scipy!")
            else:
                sys.stdout.write("This error might be caused by outdated version of scipy!\n")
            raise e
    if verbose_log:
        if log_handler:
            log_handler.info(str(results))
        else:
            sys.stdout.write(str(results) + "\n")
    best_scheme = sorted(results, key=lambda x: x["bic"])[0]
    return best_scheme


def find_greatest_common_divisor(number_list):  # euclid_algorithm
    number_list = number_list[:]
    if len(number_list) == 1:
        return number_list[0]
    elif len(number_list) == 0:
        return
    else:
        a = number_list[0]
        for i in range(len(number_list) - 1):
            a = number_list[i]
            b = number_list[i + 1]
            while b:
                a, b = b, a % b
            number_list[i + 1] = a
        return a


def reduce_list_with_gcd(number_list):
    if len(number_list) == 1:
        return [1] if number_list[0] != 0 else number_list
    elif len(number_list) == 0:
        return []
    else:
        gcd_num = find_greatest_common_divisor(number_list)
        return [int(raw_number / gcd_num) for raw_number in number_list]


