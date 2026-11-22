import argparse
import pandas as pd
import numpy as np
from scipy.stats import shapiro
import os
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import pingouin as pg
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import DistanceMatrix, permanova, permdisp
from statsmodels.multivariate.manova import MANOVA


#FUNKCJE GLOBALNE

#weryfikacja założenia normalności rozkładu
def check_normality(df, group1_cols, group2_cols, input_file):
    normality_results = []

    for i, row in df.iterrows():
        gene = row['gen']

        group1 = row[group1_cols].values.astype(float)
        group2 = row[group2_cols].values.astype(float)

        W_1, p_1 = shapiro(group1)
        W_2, p_2 = shapiro(group2)

        normality_results.append({
            'gen': gene,
            'W_1': W_1,
            'p_1': p_1,
            'W_2': W_2,
            'p_2': p_2,
            'grupa1_normal': p_1 > 0.05,
            'grupa2_normal': p_2 > 0.05
        })

    normality_df = pd.DataFrame(normality_results)

    pct_group1 = normality_df['grupa1_normal'].mean()
    pct_group2 = normality_df['grupa2_normal'].mean()

    summary_row = {
        'gen': '__podsumowanie__',
        'W_1': '',
        'p_1': '',
        'W_2': '',
        'p_2': '',
        'grupa1_normal': pct_group1,
        'grupa2_normal': pct_group2
    }

    normality_df = pd.concat(
        [normality_df, pd.DataFrame([summary_row])],
        ignore_index=True
    )

    base_name = os.path.splitext(input_file)[0]
    normality_output = base_name + "_normalnosc.csv"

    normality_df.to_csv(normality_output, index=False)

    print("Zapisano wyniki normalności do:", normality_output)
    print("Procent genów normalnych w grupie 1:", pct_group1 * 100)
    print("Procent genów normalnych w grupie 2:", pct_group2 * 100)

    return pct_group1, pct_group2

#interpretacja p-value
def interpret_p(p):
    if p < 0.05:
        return "istotne"
    elif p < 0.1:
        return "trend"
    else:
        return "nieistotne"

#kierunek log2fc
def direction(fc):
    if fc > 0:
        return "wzrost"
    elif fc < 0:
        return "spadek"
    else:
        return "brak zmiany"

#wielkość log2fc
def log2fc_size(fc):
    fc = abs(fc)
    if fc < 0.32:
        return "bardzo mały"
    elif fc < 0.58:
        return "mały"
    elif fc < 1:
        return "umiarkowany"
    else:
        return "duży"

#wielkość r_rank
def rrank_size(r):
    r = abs(r)
    if r < 0.1:
        return "bardzo mały"
    elif r < 0.3:
        return "mały"
    elif r < 0.5:
        return "umiarkowany"
    else:
        return "duży"

#interpretacja końcowa testu Manna Whitney'a
def final_interpretation_mannwhitney(row):
    p = row["p_fdr"]
    fc = row["log2FC"]
    r = row["r_rank"]
    size = row["r_rank_wielkość"]

    if p < 0.05:
          if fc > 0:
              return f"Istotny {size} wzrost"
          else:
              return f"Istotny {size} spadek"

    elif p < 0.1:
        if fc > 0:
            return f"Trend {size} wzrostu"
        else:
            return f"Trend {size} spadku"

    elif abs(r) >= 0.3 and abs(fc) >= 0.58:
        if fc > 0:
            return "Potencjalny wzrost"
        else:
            return "Potencjalny spadek"

    else:
        return "Bez zmian"

#wielkość hedges' g
def hedges_size(g):
    g = abs(g)
    if g < 0.2:
        return "bardzo mały"
    elif g < 0.5:
        return "mały"
    elif g < 0.8:
        return "średni"
    else:
        return "duży"

#wielkość eta-kwadrat
def eta_size(e):
    if e < 0.01:
        return "bardzo mały"
    elif e < 0.06:
        return "mały"
    elif e < 0.14:
        return "średni"
    else:
        return "duży"

#interpretacja końcowa test t-welch
def final_interpretation_welch(row):
    p = row["p_fdr"]
    fc = row["log2FC"]
    g = row["hedges_g"]
    eta = row["eta_kwadrat"]

    size_g = hedges_size(g)

    if p < 0.05:
        if fc > 0:
            return f"Istotny {size_g} wzrost"
        else:
            return f"Istotny {size_g} spadek"

    elif p < 0.1:
        if fc > 0:
            return f"Trend {size_g} wzrostu"
        else:
            return f"Trend {size_g} spadku"

    elif (abs(g) >= 0.5 or eta >= 0.06) and abs(fc) >= 0.58:
        if fc > 0:
            return "Potencjalny wzrost"
        else:
            return "Potencjalny spadek"

    else:
        return "Bez zmian"

#interpretacja końcowa PERMANOVA
def interpret_permanova(p_perm, p_disp):
    if p_perm < 0.05 and p_disp >= 0.05:
        return "Istotna separacja grup"

    if p_perm < 0.05 and p_disp < 0.05:
        return "Istotna separacja grup, naruszone założenie jednorodności dyspersji"

    if p_perm >= 0.05 and p_disp < 0.05:
        return "Brak separacji grup, naruszone założenie jednorodności dyspersji"

    return "Brak istotnych różnic między grupami"

#interpretacja p-value PEMRANOVA
def perm_p_interp(p):
    return "istotne" if p < 0.05 else "nieistotne"

#interpretacja p-value PERMDISP
def disp_p_interp(p):
    return "naruszone" if p < 0.05 else "spełnione"

#interpretacja p-value MANOVA
def manova_p_interp(p):
    return "istotne" if p < 0.05 else "nieistotne"

#interpretacja p-value BOX M
def box_p_interp(p):
    return "naruszone" if p < 0.05 else "spełnione"

#interpretacja końcowa MANOVA
def interpret_manova(p_m, p_box):
    if p_m < 0.05 and p_box >= 0.05:
        return "Istotne różnice między grupami"

    if p_m < 0.05 and p_box < 0.05:
        return "Istotne różnice, naruszone założenie równości macierzy kowariancji"

    if p_m >= 0.05 and p_box < 0.05:
        return "Brak różnic, naruszone założenie równości macierzy kowariancji"

    return "Brak istotnych różnic między grupami"


#GŁÓWNA FUNKCJA PROGRAMU

def run_pipeline(input_file, output_file, group1_size, transform, mode):

    # Wczytanie danych z CSV
    df = pd.read_csv(
        input_file,
        sep=';',
        decimal=',',
        skip_blank_lines=True
    )

    # Usuwanie pustych wierszy i kolumn
    df = df.dropna(axis=1, how='all')
    df = df.dropna(how="all")

    print("Rozmiar tabeli:", df.shape)
    print(df.head())

    all_gsm = [col for col in df.columns if col.startswith('GSM')]

    if group1_size > len(all_gsm):
        raise ValueError("Rozmiar grupy 1 większy niż liczba próbek")
    if len(all_gsm) == 0:
        raise ValueError("Brak kolumn GSM")

    # Transformacje
    if transform == "log2":
        print("Przeprowadzono transformację log2")
        kolumny_gsm = df.columns.str.startswith('GSM')
        df.loc[:, kolumny_gsm] = np.log2(df.loc[:, kolumny_gsm] + 1)
    else:
        print("Nie przeprowadzono transformacji log2")

    group1_cols = all_gsm[:group1_size]
    group2_cols = all_gsm[group1_size:]
    group1_size = len(group1_cols)
    group2_size = len(group2_cols)

    print("Grupa 1:", group1_size)
    print("Grupa 2:", group2_size)


    # Tryb analizy
    if mode == "single_gene":

        #wybór testu
        test = None
        if (group1_size < 10 or group2_size < 10):
            test = "mannwhitney"
        else:
            pct_group1, pct_group2 = check_normality(df, group1_cols, group2_cols, input_file)

            if (pct_group1 >= 0.8 and pct_group2 >= 0.8):
                test = "welch"
            else:
                test = "mannwhitney"

        print("Wybrany test:", test)

        if test == "mannwhitney":
            results = []

            for i, row in df.iterrows():
                gene = row['gen']
                probe = row['sonda']

                group1 = row[group1_cols].values.astype(float)
                group2 = row[group2_cols].values.astype(float)

                U, p = mannwhitneyu(group1, group2, alternative="two-sided")

                # współczynnik korelacji rangowo-biseryjnej
                r_rank = (2 * U) / (group1_size * group2_size) - 1

                mean_1 = np.mean(group1)
                mean_2 = np.mean(group2)

                median_1 = np.median(group1)
                median_2 = np.median(group2)

                # log2FC (różnica median w skali log2)
                log2fc = median_1 - median_2

                results.append({
                    "sonda": probe,
                    "gen": gene,
                    "U": U,
                    "p": p,
                    "średnia_1": mean_1,
                    "średnia_2": mean_2,
                    "mediana_1": median_1,
                    "mediana_2": median_2,
                    "n_1": group1_size,
                    "n_2": group2_size,
                    "r_rank": r_rank,
                    "log2FC": log2fc
                })

            results_df = pd.DataFrame(results)
            results_df["test"] = "Mann-Whitney"

            # korekty p-value
            results_df["p_fdr"] = multipletests(results_df["p"], method="fdr_bh")[1]
            results_df["p_bonf"] = multipletests(results_df["p"], method="bonf")[1]

            results_df = results_df.sort_values("p")

            # interpretacja p-value (FDR)
            results_df["p_fdr_interp"] = results_df["p_fdr"].apply(interpret_p)

            # kierunek log2FC
            results_df["kierunek"] = results_df["log2FC"].apply(direction)

            # wielkość log2FC
            results_df["log2fc_wielkość"] = results_df["log2FC"].apply(log2fc_size)

            # wielkość r_rank
            results_df["r_rank_wielkość"] = results_df["r_rank"].apply(rrank_size)

            # interpretacja końcowa
            results_df["końcowa_interpretacja"] = results_df.apply(final_interpretation_mannwhitney, axis=1)

        elif test == "welch":

            results = []

            for i, row in df.iterrows():
                gene = row['gen']
                probe = row['sonda']

                group1 = row[group1_cols].values.astype(float)
                group2 = row[group2_cols].values.astype(float)

                t_stat, p = stats.ttest_ind(group1, group2, equal_var=False, alternative="two-sided")

                # hedges' g
                hedges = pg.compute_effsize(group1, group2, eftype='hedges')

                mean_1 = np.mean(group1)
                mean_2 = np.mean(group2)

                # eta-kwadrat (η²)

                var1 = np.var(group1, ddof=1)
                var2 = np.var(group2, ddof=1)

                df_welch = ((var1/group1_size + var2/group2_size)**2) / \
                          ((var1/group1_size)**2/(group1_size-1) + (var2/group2_size)**2/(group2_size-1))

                eta_sq = (t_stat**2) / (t_stat**2 + df_welch)

                # log2FC (różnica średnich w skali log2)
                log2fc = mean_1 - mean_2

                results.append({
                    "sonda": probe,
                    "gen": gene,
                    "t_stat": t_stat,
                    "p": p,
                    "średnia_1": mean_1,
                    "średnia_2": mean_2,
                    "std_1": np.std(group1, ddof=1),
                    "std_2": np.std(group2, ddof=1),
                    "n_1": group1_size,
                    "n_2": group2_size,
                    "hedges_g": hedges,
                    "eta_kwadrat": eta_sq,
                    "log2FC": log2fc
                })

            results_df = pd.DataFrame(results)
            results_df["test"] = "t-Welch"

            # Korekty p-value
            results_df["p_fdr"] = multipletests(results_df["p"], method="fdr_bh")[1]
            results_df["p_bonf"] = multipletests(results_df["p"], method="bonf")[1]

            results_df = results_df.sort_values("p")

            # interpretacja p-value (FDR)
            results_df["p_fdr_interp"] = results_df["p_fdr"].apply(interpret_p)

            # kierunek log2FC
            results_df["kierunek"] = results_df["log2FC"].apply(direction)

            # wielkość log2FC
            results_df["log2fc_wielkość"] = results_df["log2FC"].apply(log2fc_size)

            # wielkość hedges'g
            results_df["hedges_g_wielkość"] = results_df["hedges_g"].apply(hedges_size)

            # wielkość eta_squared
            results_df["eta_wielkość"] = results_df["eta_kwadrat"].apply(eta_size)

            # interpretacja końcowa
            results_df["końcowa_interpretacja"] = results_df.apply(final_interpretation_welch, axis=1)

        else:
            raise ValueError(f"Nieznany test: {test}")

    elif mode == "signature":

        #wybór testu
        test = None

        p = df.shape[0]

        n_min = min(group1_size, group2_size)

        if n_min < 15:
            test = "permanova"
        elif n_min < 5 * p:
            test = "permanova"
        else:
            pct_group1, pct_group2 = check_normality(df, group1_cols, group2_cols, input_file)

            if (pct_group1 >= 0.8 and pct_group2 >= 0.8):
                test = "manova"
            else:
                test = "permanova"

        print("Wybrany test:", test)

        if test == "permanova":

            X = df[group1_cols + group2_cols].copy()

            X = X.T

            groups = (["group1"] * group1_size + ["group2"] * group2_size)

            dist = pdist(X.values, metric="euclidean")
            dm = DistanceMatrix(
                squareform(dist),
                ids=[str(i) for i in range(len(groups))]
            )

            # PERMANOVA
            permanova_res = permanova(dm, grouping=groups, permutations=999)

            # PERMDISP
            permdisp_res = permdisp(dm, grouping=groups, permutations=999)

            n = permanova_res["sample size"]
            g = permanova_res["number of groups"]
            pseudo_f = permanova_res["test statistic"]

            results = [{
                "test": "permanova",
                "n_1": group1_size,
                "n_2": group2_size,
                "l_zmiennych": p,
                "pseudo_F": pseudo_f,
                "permanova_p": permanova_res["p-value"],
                "permdisp_F": permdisp_res["test statistic"],
                "permdisp_p": permdisp_res["p-value"],
                "p_perm_interp": perm_p_interp(permanova_res["p-value"]),
                "p_disp_interp": disp_p_interp(permdisp_res["p-value"]),
                "końcowa_interpretacja": interpret_permanova(permanova_res["p-value"], permdisp_res["p-value"])
            }]

            results_df = pd.DataFrame(results)


        elif test == "manova":

            X = df[group1_cols + group2_cols].copy()

            X = X.T
            X.index = range(X.shape[0])

            groups = (["group1"] * group1_size + ["group2"] * group2_size)

            X.columns = ["gene_" + str(i) for i in X.columns]
            X["group"] = groups

            dependent_vars = X.columns[:-1]
            formula = " + ".join(dependent_vars) + " ~ group"

            # MANOVA
            mv = MANOVA.from_formula(formula, data=X)
            mv_res = mv.mv_test()

            stat_table = mv_res.results["group"]["stat"]

            # BOX M test
            box_res = pg.box_m(
                data=X,
                dvs=dependent_vars,
                group="group"
            )

            p_m = stat_table.loc["Pillai's trace", "Pr > F"]
            p_box = box_res["pval"].values[0]

            results = [{
                "test": "manova",
                "manova_typ": "pillai",
                "n_1": group1_size,
                "n_2": group2_size,
                "l_zmiennych": len(dependent_vars),
                "pillai_trace": stat_table.loc["Pillai's trace", "Value"],
                "manova_F": stat_table.loc["Pillai's trace", "F Value"],
                "manova_p": p_m,
                "box_m_stat": box_res["Chi2"].values[0],
                "box_m_df": box_res["df"].values[0],
                "box_m_p": p_box,
                "manova_p_interp": manova_p_interp(p_m),
                "box_p_interp": box_p_interp(p_box),
                "końcowa_interpretacja": interpret_manova(p_m, p_box)
            }]

            results_df = pd.DataFrame(results)

        else:
            raise ValueError(f"Nieznany test: {test}")


    # Zapis do pliku CSV
    results_df.to_csv(output_file, index=False)
    print(f"Wyniki zapisane w {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Program do automatycznej analizy danych ekspresji genów "
        ),
        epilog=(
            "Przykład użycia:\n"
            "python program.py -i dane.csv -o wyniki.csv "
            "--group1_size 5 --transform log2 --mode single_gene"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-i", "--input", required=True, help="nazwa pliku CSV wejściowego")
    parser.add_argument("-o", "--output", required=True, help="nazwa pliku wynikowego")
    parser.add_argument("--group1_size", type=int, default=1, help="liczba próbek w grupie 1")
    parser.add_argument("--transform", choices=["none", "log2"], default="none", help="tranformacja danych przed analizą")
    parser.add_argument("--mode", choices=["single_gene", "signature"], default="single_gene",help=(
        "tryb analizy:\n"
        "single_gene  - analiza pojedynczych genów\n"
        "signature    - analiza sygnatur genowych"
    ))


    args = parser.parse_args()
    run_pipeline(args.input, args.output,args.group1_size, args.transform, args.mode)

if __name__ == "__main__":
    main()
