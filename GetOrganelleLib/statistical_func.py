# try:
#     from scipy import stats, inf, log
# except ImportError:
#     class stats:
#         class norm:
#             def logpdf(foo1, foo2, foo3):
#                 raise ImportError("Failed in 'from scipy import stats, inf, log'!")
#     inf = float("inf")
from math import log, inf, sqrt, pi
from itertools import permutations
from copy import deepcopy
# add try except so that when statistical_func.py is called by other scripts, it will not prompt error immediately
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
import sys
try:
    from gekko import GEKKO
except ImportError:
    def GEKKO(remote):
        raise ImportError("Failed in 'from gekko import GEKKO'!")


def bic(loglike, len_param, len_data):
    return log(len_data) * len_param - 2 * loglike


def aic(loglike, len_param):
    return 2 * len_param - 2 * loglike


def weighted_mean_and_std(values, weights):
    mean = np.average(values, weights=weights)
    std = np.average((values-mean)**2, weights=weights)**0.5
    return mean, std


def norm_logpdf(numpy_array, mu, sigma):
    u = (numpy_array-mu)/abs(sigma)
    return log(1/(sqrt(2*pi)*abs(sigma)))-u*u/2


def weighted_clustering_with_em_aic(data_array,
                                    data_weights=None,
                                    minimum_cluster=1,
                                    maximum_cluster=5,
                                    min_sigma_factor=1E-5,
                                    cluster_limited=None,
                                    cluster_bans=None,
                                    log_handler=None,
                                    verbose_log=False,
                                    random_obj=None):
    """
    The current implementation is using a categorical distribution,
        with assignment of data exclusively to specific components.

    :param data_array:
    :param data_weights:
    :param minimum_cluster:
    :param maximum_cluster:
    :param min_sigma_factor:
    :param cluster_limited: {dat_id1: {0, 1}, dat_id2: {0}, dat_id3: {0} ...}
           cluster_limited has priority over cluster_bans if conflicted
    :param cluster_bans: {dat_a: {0, 1}, dat_b: {0}, dat_c: {0} ...}
    :param log_handler:
    :param verbose_log:
    :param random_obj: to control the random process using a universal seed
    :return:
    """
    # import time
    # time0 = time.time()
    # pdf_time = [0.]
    if random_obj is None:
        import random as random_obj

    min_sigma = min_sigma_factor * np.average(data_array, weights=data_weights)

    def model_loglike(dat_arr, dat_w, lbs, parameters):
        total_loglike = 0
        for go_to_cl, pr in enumerate(parameters):
            points = dat_arr[lbs == go_to_cl]
            weights = dat_w[lbs == go_to_cl]
            if len(points):
                # total_loglike += sum(norm_logpdf(points, pr["mu"], pr["sigma"]) * weights)
                # total_loglike += sum(norm_logpdf(points, pr["mu"], pr["sigma"]) * weights + log(pr["percent"]))
                # total_loglike += sum(norm_logpdf(points, pr["mu"], pr["sigma"]) * weights + log(pr["percent"]))
                total_loglike += sum((norm_logpdf(points, pr["mu"], pr["sigma"]) + log(pr["percent"])) * weights)
        return total_loglike

    def revise_labels_according_to_constraints(_raw_lbs, _fixed, loglike_table=None):
        _new_labels = list(_raw_lbs)
        for _dat_id in range(int(data_len)):
            if _dat_id in _fixed:
                if _raw_lbs[_dat_id] in _fixed[_dat_id]:
                    # _new_labels[_dat_id] = _raw_lbs[_dat_id]
                    pass
                else:
                    if loglike_table is None:
                        _new_labels[_dat_id] = sorted(_fixed[_dat_id])[0]
                    else:
                        val_increasing_order = list(loglike_table[:, _dat_id].argsort())
                        for order_id in range(len(val_increasing_order) - 2, -1, -1):  # search from the sub-optimum
                            this_label = val_increasing_order.index(order_id)
                            if this_label in _fixed[_dat_id]:
                                _new_labels[_dat_id] = this_label
                                # v1.8.0-pre6 fix a bug
                                break
            else:
                # _new_labels[_dat_id] = _raw_lbs[_dat_id]
                pass
        return np.array(_new_labels)

    def assign_cluster_labels(dat_arr, dat_w, parameters, lb_fixed):
        # assign every data point to its most likely cluster
        if len(parameters) == 1:
            return np.array([0] * int(data_len))
        # elif len(in_params) == len(dat_arr):
        #     return np.array(range(len(in_params)))
        else:
            # the parameter set of the first cluster
            # timex = time.time()
            # loglike_res = norm_logpdf(dat_arr, parameters[0]["mu"], parameters[0]["sigma"]) * dat_w
            # loglike_res = norm_logpdf(dat_arr, parameters[0]["mu"], parameters[0]["sigma"]) * dat_w + \
            #               log(parameters[0]["percent"])
            loglike_res = (norm_logpdf(dat_arr, parameters[0]["mu"], parameters[0]["sigma"]) +
                           log(parameters[0]["percent"])) * dat_w
            # the parameter set of the rest cluster
            for pr in parameters[1:]:
                loglike_res = np.vstack(
                    # (loglike_res, norm_logpdf(dat_arr, pr["mu"], pr["sigma"]) * dat_w))
                    # (loglike_res, norm_logpdf(dat_arr, pr["mu"], pr["sigma"]) * dat_w + log(pr["percent"])))
                    (loglike_res, (norm_logpdf(dat_arr, pr["mu"], pr["sigma"]) + log(pr["percent"])) * dat_w))
            # print(loglike_res)
            # pdf_time[0] += time.time() - timex
            # assign labels
            new_labels = loglike_res.argmax(axis=0)
            if lb_fixed:
                # intermediate_labels = []
                # for here_dat_id in range(int(data_len)):
                #     if here_dat_id in lb_fixed:
                #         if new_labels[here_dat_id] in lb_fixed[here_dat_id]:
                #             intermediate_labels.append(new_labels[here_dat_id])
                #         else:
                #             intermediate_labels.append(sorted(lb_fixed[here_dat_id])[0])
                #     else:
                #         intermediate_labels.append(new_labels[here_dat_id])
                # np.array(intermediate_labels)
                new_labels = revise_labels_according_to_constraints(new_labels, lb_fixed, loglike_res)
                # new_labels = np.array([
                # sorted(cluster_limited[dat_item])[0]
                # if new_labels[here_dat_id] not in cluster_limited[dat_item] else new_labels[here_dat_id]
                # if dat_item in cluster_limited else
                # new_labels[here_dat_id]
                # for here_dat_id, dat_item in enumerate(data_array)])
            #     limited_values = set(dat_arr[list(lb_fixed)])
            # else:
            #     limited_values = set()
            if len(set(new_labels)) > len(parameters):
                raise ValueError("Assigning failed!")
            label_counts = {lb: 0 for lb in range(len(parameters))}
            for ct_lb in new_labels:
                label_counts[ct_lb] += 1
            # 2023-01-15 added
            empty_lbs = []
            # filled_lbs = []
            for label_id in range(len(parameters)):
                if label_counts[label_id] == 0:
                    empty_lbs.append(label_id)
                # else:
                #     filled_lbs.append(label_id)
            if empty_lbs:
                # find the new option that reduce the likelihood least
                # if len(empty_lbs) == 1:
                #     new_lb = empty_lbs[0]
                #     loglike_reduced = np.max(loglike_res, axis=0) - loglike_res[new_lb]
                #     orders_to_d = {order_id: _d_id for _d_id, order_id in enumerate(loglike_reduced.argsort())}
                #     for try_order_id in range(len(dat_arr)):
                #         try_data_id = orders_to_d[try_order_id]
                #         # if old label is not fixed, or the new label is within the constraint
                #         if try_data_id not in lb_fixed or new_lb in lb_fixed[try_data_id]:
                #             new_labels[try_data_id] = new_lb
                #             break
                #     else:
                #         raise ValueError("Assigning failed!")
                # else:
                tmp_info = {}
                to_change = {}
                lb_counts = deepcopy(label_counts)
                for new_lb in empty_lbs:
                    tmp_info[new_lb] = {}
                    loglike_reduced = np.max(loglike_res, axis=0) - loglike_res[new_lb]
                    tmp_info[new_lb]["reduced"] = loglike_reduced
                    orders_to_d = {order_id: _d_id for _d_id, order_id in enumerate(loglike_reduced.argsort())}
                    tmp_info[new_lb]["orders"] = orders_to_d
                    for try_order_id in range(len(dat_arr)):
                        try_data_id = orders_to_d[try_order_id]
                        # if old label is not fixed, or the new label is within the constraint
                        if (try_data_id not in lb_fixed or new_lb in lb_fixed[try_data_id]) and \
                                lb_counts[new_labels[try_data_id]] > 1:
                            to_change[try_data_id] = new_lb
                            lb_counts[new_lb] += 1
                            lb_counts[new_labels[try_data_id]] -= 1
                            break
                    else:
                        raise ValueError("Assigning failed!")
                if len(to_change) == len(empty_lbs):  # each empty lbs has a unique data to fill
                    for data_id, new_lb in to_change.items():
                        new_labels[data_id] = new_lb
                else:
                    # slow but easy to coding way
                    replace_res = {}
                    for new_lb_order in permutations(empty_lbs):
                        replace_res[new_lb_order] = {"loglike": 0.,
                                                     "change": {},
                                                     "failed": False,
                                                     "counts": deepcopy(label_counts)}
                        for new_lb in new_lb_order:
                            for try_order_id in range(len(dat_arr)):
                                try_data_id = tmp_info[new_lb]["orders"][try_order_id]
                                # if the data id was not used "by other new lbs", AND
                                # if the donor has more than two occurrences
                                # if data id is not fixed, or the new label is within the constraint
                                if try_data_id not in replace_res[new_lb_order]["change"] and \
                                        replace_res[new_lb_order]["counts"][new_labels[try_data_id]] > 1 and \
                                        (try_data_id not in lb_fixed or new_lb in lb_fixed[try_data_id]):
                                    replace_res[new_lb_order]["change"][try_data_id] = new_lb
                                    replace_res[new_lb_order]["loglike"] -= tmp_info[new_lb]["reduced"][try_data_id]
                                    replace_res[new_lb_order]["counts"][new_lb] += 1
                                    replace_res[new_lb_order]["counts"][new_labels[try_data_id]] -= 1
                                    break
                            else:
                                replace_res[new_lb_order]["failed"] = True
                                replace_res[new_lb_order]["loglike"] = -inf
                    best_order, best_info = sorted(replace_res.items(), key=lambda x: -x[1]["loglike"])[0]
                    if best_info["failed"]:
                        raise ValueError("Assigning failed!")
                    else:
                        for data_id, new_lb in best_info["change"].items():
                            new_labels[data_id] = new_lb
            # if there is an empty cluster,
            # and if there is another non-empty cluster with two ends not in the fixed (lb_fixed),
            # then move one of the end (min or max) from that non-empty cluster to the empty cluster
            # for empty_lb in label_counts:
            #     if label_counts[empty_lb] == 0:
            #         non_empty_lbs = {ne_lb: [min, max] for ne_lb in label_counts if label_counts[ne_lb] > 1}
            #         for af_lb in sorted(non_empty_lbs):
            #             these_points = dat_arr[new_labels == af_lb]
            #             if max(these_points) in limited_values:
            #                 non_empty_lbs[af_lb].remove(max)
            #             if min(these_points) in limited_values:
            #                 non_empty_lbs[af_lb].remove(min)
            #             if not non_empty_lbs[af_lb]:
            #                 del non_empty_lbs[af_lb]
            #         if non_empty_lbs:
            #             chose_lb = random_obj.choice(list(non_empty_lbs))
            #             chose_points = dat_arr[new_labels == chose_lb]
            #             # random.choice([min, max]), then use the resulting function to pick the point
            #             data_point = random_obj.choice(non_empty_lbs[chose_lb])(chose_points)
            #             # random.choice(np.array([0])) triggers: IndexError: Cannot choose from an empty sequence
            #             transfer_index = random_obj.choice(list(np.where(dat_arr == data_point)[0]))
            #             new_labels[transfer_index] = empty_lb
            #             label_counts[chose_lb] -= 1
            #             # 2022-12-18 fix a long-lasting issue
            #             label_counts[empty_lb] += 1
            return new_labels

    def updating_parameter(dat_arr, dat_w, lbs, in_params):
        new_params = deepcopy(in_params)
        for go_to_cl, pr in enumerate(new_params):
            these_points = dat_arr[lbs == go_to_cl]
            these_weights = dat_w[lbs == go_to_cl]
            if len(these_points) > 1:
                this_mean, this_std = weighted_mean_and_std(these_points, these_weights)
                pr["mu"] = this_mean
                pr["sigma"] = max(this_std, min_sigma)
                pr["percent"] = sum(these_weights)  # / data_len
            elif len(these_points) == 1:
                pr["sigma"] = max(dat_arr.std() / data_len, min_sigma)
                # 2023-01-15
                # pr["mu"] = np.average(these_points, weights=these_weights) + pr["sigma"] * (2 * random_obj.random() - 1)
                pr["mu"] = these_points[0]
                pr["percent"] = sum(these_weights)  # / data_len
            else:
                # exclude
                pr["mu"] = max(dat_arr) * 1E4
                pr["sigma"] = min_sigma
                pr["percent"] = 1E-10
        return new_params

    data_array = np.array(data_array)
    data_len = float(len(data_array))
    if data_weights is None or not len(data_weights):
        data_weights = np.array([1. for foo in range(int(data_len))])
    else:
        assert len(data_weights) == data_len
        average_weights = float(sum(data_weights)) / data_len
        # normalized
        data_weights = np.array([raw_w / average_weights for raw_w in data_weights])

    results = []
    # adjust the min and max number of clusters according to constraints
    if cluster_limited is None:
    #     freedom_dat_item = int(data_len)
        cluster_limited = {}
    # else:
    #     cls = set()
    #     for sub_cls in cluster_limited.values():
    #         cls |= sub_cls
    #     freedom_dat_item = max(0, int(data_len) - max(0, len(cluster_limited) - len(cls)))
    # min_choices = 0
    # if cluster_bans is None:
    #     cluster_bans = {}
    # else:
    #     for ban_dat_id, sub_cls in cluster_bans.items():
    #         if ban_dat_id not in cluster_limited:  # cluster_limited has priority over cluster_bans
    #             min_choices = max(len(sub_cls) + 1, min_choices)
    # # assert min_choices < freedom_dat_item, "unrealistic constraints: \ncluster_limited: " + \
    # #                                        str(cluster_limited) + "\ncluster_bans: " + \
    # #                                        str(cluster_bans)
    # minimum_cluster = min(freedom_dat_item, max(minimum_cluster, min_choices))
    # maximum_cluster = min(freedom_dat_item, max(maximum_cluster, min_choices))
    # 2023-01-15
    # print(minimum_cluster, maximum_cluster)
    minimum_cluster = max(minimum_cluster, len(set([tuple(_cl) for _cl in cluster_limited.values() if len(_cl) == 1])))
    maximum_cluster = min(maximum_cluster, len(data_array))
    # print(minimum_cluster, maximum_cluster)

    # timey = time.time()
    # round_times = []
    # iteratively try the num of clusters
    for total_cluster_num in range(minimum_cluster, maximum_cluster + 1):
        cluster_num_failure = False
        if log_handler and verbose_log:
            log_handler.info("assessing %i clusters" % total_cluster_num)
        if cluster_limited:
            this_limit = deepcopy(cluster_limited)
            for dat_id in cluster_bans:
                if dat_id not in this_limit:
                    this_limit[dat_id] = set()
                    for potential_lb_id in range(total_cluster_num):
                        if potential_lb_id not in cluster_bans[dat_id]:
                            this_limit[dat_id].add(potential_lb_id)
        else:
            this_limit = {}
        # initialization
        # # labels = np_rd_obj.choice(total_cluster_num, int(data_len))
        # # TODO: each cluster has to have at least one occurrence, which will be complicated combined with this_limit
        # min_occurrences = list(range(total_cluster_num))
        # random_obj.shuffle(min_occurrences)
        # labels = np.array(min_occurrences +
        #                   random_obj.choices(range(total_cluster_num), k=int(data_len) - total_cluster_num))
        # labels = revise_labels_according_to_constraints(labels, this_limit)
        # using decreasing order rather than random @2023-01-18, to create more reasonable initial parameter set
        extra_for_first_lb = int(data_len) % total_cluster_num
        each_len = (int(data_len) - extra_for_first_lb) // total_cluster_num
        cov_decreasing_order = np.flip(data_array.argsort())
        labels = np.zeros(int(data_len), dtype=np.int8)
        for go_l, go_d in enumerate(range(extra_for_first_lb, int(data_len), each_len)):
            labels[cov_decreasing_order[go_d: go_d + each_len]] = go_l
        # initialize the parameters
        norm_parameters = updating_parameter(data_array, data_weights, labels,
                                             # [{"mu": 0, "sigma": 1}
                                             [{"mu": 0, "sigma": 1, "percent": total_cluster_num/data_len}
                                              for foo in range(total_cluster_num)])
        if log_handler and verbose_log:
            log_handler.info("    initial labels: " + str(list(labels)))
            log_handler.info("    initial params: " + str(norm_parameters))
        loglike_shift = inf
        prev_loglike = -inf
        epsilon = 0.001
        count_iterations = 0
        best_loglike = -inf
        best_parameter = norm_parameters
        count_best = 1
        while loglike_shift > epsilon and count_iterations < 500 and count_best < 50:
            count_iterations += 1
            # expectation
            try:
                labels = assign_cluster_labels(data_array, data_weights, norm_parameters, this_limit)
            except ValueError as e:
                if str(e) == "Assigning failed!":
                    if verbose_log and log_handler:
                        log_handler.info("    assigning failed for %i clusters" % total_cluster_num)
                    cluster_num_failure = True
                    break
                else:
                    raise e
            # maximization
            updated_parameters = updating_parameter(data_array, data_weights, labels, deepcopy(norm_parameters))
            this_loglike = model_loglike(data_array, data_weights, labels, updated_parameters)
            if log_handler and verbose_log:
                log_handler.info("    iter_%i labels: " % count_iterations + str(list(labels)))
                log_handler.info("    iter_%i params: " % count_iterations + str(updated_parameters))
                log_handler.info("    iter_%i loglike: " % count_iterations + str(this_loglike))
            loglike_shift = abs((this_loglike - prev_loglike) / this_loglike)
            # update
            prev_loglike = this_loglike
            norm_parameters = updated_parameters
            if this_loglike > best_loglike:
                best_parameter = updated_parameters
                best_loglike = this_loglike
                count_best = 1
            else:
                count_best += 1
        if cluster_num_failure:
            break
        # 2023-01-15 replace: labels = assign_cluster_labels(data_array, data_weights, best_parameter, None)
        labels = assign_cluster_labels(data_array, data_weights, best_parameter, this_limit)
        results.append({"loglike": best_loglike, "iterates": count_iterations, "cluster_num": total_cluster_num,
                        "parameters": best_parameter, "labels": labels,
                        "aic": aic(best_loglike, 2 * total_cluster_num),
                        "bic": bic(best_loglike, 2 * total_cluster_num, data_len)})
        # except TypeError as e:
        #     if log_handler:
        #         log_handler.error("This error might be caused by outdated version of scipy!")
        #     else:
        #         sys.stdout.write("This error might be caused by outdated version of scipy!\n")
        #     raise e
        # round_times.append(time.time() - timey)
        # timey = time.time()
    if verbose_log:
        if log_handler:
            log_handler.info(str(results))
        else:
            sys.stdout.write(str(results) + "\n")
    # TODO: if all clustering failed
    if not results:
        raise ValueError("Solution Not Found!")
    best_scheme = sorted(results, key=lambda x: x["bic"])[0]
    # print(time.time() - time0, round_times, pdf_time)
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


