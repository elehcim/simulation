import pandas as pd
from simulation.simdata import get_pickle_data, get_last_d, get_first_d, last_d2df
from functools import lru_cache
from simulation.util import setup_logger

logger = setup_logger('data_handler', logger_level='INFO')

def rtre_apo_idx_2(df_all):
    """Pick the first occurence of the tidal radius over effective radius criterion"""
    f = df_all.rt10/df_all.r_eff3d
    return f.lt(1).idxmax()


def cut_out_rt_criterion(df_big):

    # This contains the index of the first occurrence. If the criterion is never satisfied it is zero.
    ci = df_big.groupby(['name', 'pericenter'], sort=False).apply(rtre_apo_idx_2)

    # I use this dict format again with those keys because I have already some useful functions to transform them
    d = dict()
    for ((full_name, peri), group) in df_big.groupby(['name', 'pericenter'], sort=False):
        name = full_name[:2]
        stop_idx = ci[full_name][peri]
        # print(full_name, peri, stop_idx)
        key_name = f'{name}p{peri}'
        if stop_idx == 0:
            d[key_name] = group
        else:
            d[key_name] = group.iloc[:stop_idx]
    return d


class DataHandler:
    @classmethod
    def show_data_files(cls):
        """Show data files sorted by date"""
        import os
        import glob
        from simulation.simdata import _get_data_dir
        pkl_files = list(map(os.path.basename, glob.glob(os.path.join(_get_data_dir(), '*.pkl'))))
        pkl_files.sort(key=lambda x: x.split('_')[-1])
        return tuple(pkl_files)

    def __init__(self, cache_file=None):
        """
        Quickly get data from various sources, using simulation.simdata

        Parameters
        ----------
        cache_file : str
        """
        if cache_file is None:
            logger.info(f"Getting most recent cache file")
            cache_file = self.show_data_files()[-1]
        self.cache_file = cache_file

    @lru_cache(1)
    def data(self):
        return get_pickle_data(self.cache_file)

    def data_last(self):
        last_d = get_last_d(self.data())
        last_df = last_d2df(last_d)
        return last_df

    def data_first(self):
        first_d = get_first_d(self.data())
        first_df = last_d2df(first_d)
        return first_df

    @lru_cache(1)
    def data_rt(self):
        d_rt = cut_out_rt_criterion(self.data_big())
        return d_rt

    def data_big_rt(self):
         return pd.concat([v for v in self.data_rt().values()], axis=0)

    def data_last_rt(self):
        last_d = get_last_d(self.data_rt())
        last_df = last_d2df(last_d)
        return last_df


    def data_big(self):
        db = pd.concat([v for v in self.data().values()], axis=0, sort=True)
        # Uniques are returned in order of appearance.
        # Even if pd.unique does NOT sort for us, they appear in sorted order, so it's ok.
        # concat inspired from here:https://stackoverflow.com/a/35850749
        sim_label = db['name'].str.slice(stop=2).str.cat(db['pericenter'].astype(str), sep='p')
        # from here: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.astype.html
        sim_label_dtype = pd.api.types.CategoricalDtype(sim_label.unique(), ordered=True)
        db['sim_label'] = sim_label.astype(sim_label_dtype)
        name_dtype = pd.api.types.CategoricalDtype(db['name'].unique(), ordered=True)
        db['name'] = db['name'].astype(name_dtype)

        return db
