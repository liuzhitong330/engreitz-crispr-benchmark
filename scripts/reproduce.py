#!/usr/bin/env python3
"""Reproduce the saved example cohort, audit joins, and compare fixed rankings.

No fitting, threshold optimization, label creation, or untested negative sampling.
Official numerical R reproduction is separate: reproduce_official.R.
"""
import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, precision_recall_curve
from fetch_inputs import fetch, ROOT

SEED = 20261029
BINS = [('all', 0, np.inf), ('0–20 kb', 0, 20000),
        ('20–100 kb', 20000, 100000), ('100–2,500 kb', 100000, 2500000)]

def overlap_mask(start, end, starts, ends, half_open=False):
    """BED half-open for TSS filtering; original inclusive semantics for pair join."""
    if half_open:
        return (starts < end) & (ends > start)
    return (starts <= end) & (ends >= start)

def blocked_by_tss(pred, tss):
    """Half-open interval overlap via sorted TSS starts / cumulative maximum ends."""
    blocked = np.zeros(len(pred), dtype=bool)
    for chrom, idx in pred.groupby('chr').groups.items():
        windows = tss[tss['#chr'] == chrom].sort_values('start')
        if windows.empty:
            continue
        starts = windows.start.to_numpy()
        maxends = np.maximum.accumulate(windows.end.to_numpy())
        elements = pred.loc[idx]
        last = np.searchsorted(starts, elements.end.to_numpy(), side='left') - 1
        blocked[idx] = (last >= 0) & (maxends[np.maximum(last, 0)] > elements.start.to_numpy())
    return blocked

def build_cohort(inputs):
    expt = pd.read_csv(inputs / 'resources/example/EPCrisprBenchmark_Fulco2019_K562_GRCh38.tsv.gz', sep='\t', dtype={'ValidConnection': str})
    pred = pd.read_csv(inputs / 'resources/example/ABC_K562_Fulco2019Genes_GRCh38.tsv.gz', sep='\t')
    tss = pd.read_csv(inputs / 'resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed', sep='\t')
    config = pd.read_csv(inputs / 'resources/example/pred_config.txt', sep=r'\s+')
    assert config.iloc[0].pred_col == 'powerlaw.Score'
    assert config.iloc[0].aggregate_function == 'sum' and config.iloc[0].fill_value == 0
    assert float(config.iloc[0].alpha) == .02
    assert set(expt.CellType) == {'K562'} and set(pred.CellType) == {'K562_ID_2644'}
    mapping = pd.read_csv(inputs / 'resources/example/ABC_cell_type_mapping.txt', sep=r'\s+')
    assert mapping[['experiment','predictions']].values.tolist() == [['K562','K562_ID_2644']]
    pred['mapped_cell'] = pred.CellType.map(mapping.set_index('predictions').experiment)
    assert pred.mapped_cell.notna().all()
    assert not tss.name.duplicated().any(), 'Ambiguous gene symbol in supplied TSS universe'
    assert pred['powerlaw.Score'].notna().all()
    assert (pred.start < pred.end).all() and (expt.chromStart < expt.chromEnd).all()
    tss_by_gene = tss.set_index('name')
    blocked = blocked_by_tss(pred, tss)
    filtered = pred[~blocked]
    groups = {(cell, chrom, gene): p for (cell, chrom, gene), p in filtered.groupby(['mapped_cell', 'chr', 'TargetGene'])}
    raw_groups = {(cell, chrom, gene): p for (cell, chrom, gene), p in pred.groupby(['mapped_cell', 'chr', 'TargetGene'])}
    valid = expt[(expt.ValidConnection == 'TRUE') & expt.Regulated.notna() & expt.measuredGeneSymbol.isin(tss_by_gene.index)].copy()
    assert not valid.duplicated(['CellType', 'chrom', 'chromStart', 'chromEnd', 'measuredGeneSymbol']).any()
    rows = []
    endpoint_differences = 0
    for r in valid.itertuples(index=False):
        gene = r.measuredGeneSymbol
        annot = tss_by_gene.loc[gene]
        assert annot['#chr'] == r.chrom, 'Gene/enhancer chromosome mismatch'
        center = int((int(annot.start) + int(annot.end)) // 2)
        gene_predictions = pred[pred.TargetGene == gene]
        assert gene_predictions.TargetGeneTSS.eq(center).all(), 'Prediction/gene-universe TSS mismatch'
        key = (r.CellType, r.chrom, gene)
        matches = groups.get(key, filtered.iloc[0:0])
        mask = overlap_mask(r.chromStart, r.chromEnd, matches.start, matches.end)
        half = overlap_mask(r.chromStart, r.chromEnd, matches.start, matches.end, True)
        endpoint_differences += int((mask != half).sum())
        hits = matches[mask]
        raw = raw_groups.get(key, pred.iloc[0:0])
        raw_hits = raw[overlap_mask(r.chromStart, r.chromEnd, raw.start, raw.end)]
        score = float(hits['powerlaw.Score'].sum()) if len(hits) else 0.
        regulated = bool(r.Regulated)
        category = ('ABC high / assay down' if regulated else 'ABC high / no down call') if score >= .02 else ('ABC low / assay down' if regulated else 'ABC low / no down call')
        rows.append(dict(id=r.name, gene=gene, chrom=r.chrom, start=int(r.chromStart), end=int(r.chromEnd),
                         tss=center, original_expt_tss=r.startTSS,
                         distance_bp=max(abs((int(r.chromStart)+int(r.chromEnd))//2-center)-1, 0),
                         absolute_midpoint_distance_bp=abs((int(r.chromStart)+int(r.chromEnd))/2-center),
                         abc_score=score, matched=bool(len(hits)), overlap_count=len(hits), raw_overlap_count=len(raw_hits),
                         regulated=regulated, significant=bool(r.Significant), effect=float(r.EffectSize),
                         padj=r.pValueAdjusted, power=r.PowerAtEffectSize25, reference=r.Reference, cell_type=r.CellType,
                         missing_reason=('covered' if len(hits) else 'TSS-filtered only' if len(raw_hits) else 'absent before TSS filter'),
                         pred_elements=';'.join(sorted(hits.name.astype(str))), category=category))
    cohort = pd.DataFrame(rows)
    excluded = expt[~expt.name.isin(cohort.id)].copy()
    excluded['exclusion_reason'] = np.where(~excluded.measuredGeneSymbol.isin(tss_by_gene.index), 'gene absent from supplied universe',
                                           np.where(excluded.ValidConnection != 'TRUE', excluded.ValidConnection, 'missing Regulated label'))
    assert cohort.regulated.tolist() == (valid.Significant & (valid.EffectSize < 0)).tolist()
    assert not cohort.id.duplicated().any()
    embedded = valid.name.str.extract(r'\|([^:]+):(\d+)-(\d+):')
    historical_names = ((embedded[0].to_numpy() != valid.chrom.to_numpy()) |
                        (pd.to_numeric(embedded[1]).to_numpy() != valid.chromStart.to_numpy()) |
                        (pd.to_numeric(embedded[2]).to_numpy() != valid.chromEnd.to_numpy())).sum()
    qc = dict(experimental_rows=len(expt), prediction_rows=len(pred), prediction_rows_removed=int(blocked.sum()),
              prediction_rows_retained=len(filtered), valid_pairs=len(cohort), positives=int(cohort.regulated.sum()), genes=cohort.gene.nunique(),
              missing=int((~cohort.matched).sum()), missing_positives=int((~cohort.matched & cohort.regulated).sum()),
              missing_by_reason=cohort[~cohort.matched].missing_reason.value_counts().to_dict(),
              overlap_multiplicity={str(k): int(v) for k,v in cohort.overlap_count.value_counts().sort_index().items()},
              endpoint_overlap_differences=endpoint_differences,
              name_coordinate_disagreements=int(historical_names),
              genes_original_tss_different=int(cohort.loc[cohort.original_expt_tss != cohort.tss, 'gene'].nunique()),
              prediction_ensembl_all_missing=bool(pred.TargetGeneEnsemblID.isna().all()),
              prediction_distance_all_missing=bool(pred.DistanceToTSS.isna().all()),
              negative_significant_increases=int((~cohort.regulated & cohort.significant & (cohort.effect > 0)).sum()),
              negative_power_below_80pct=int((~cohort.regulated & (cohort.power < .8)).sum()),
              experimental_references=expt.Reference.value_counts().to_dict(),
              evaluated_references=cohort.reference.value_counts().to_dict(),
              original_label_mismatches=0, genome='GRCh38', cell_mapping='K562_ID_2644 → K562',
              abc_threshold=.02, gene_identifier_join='Exact symbols; no alias or Ensembl-ID remapping',
              distance_convention='Upstream GRanges midpoint-to-TSS intervening bases; 1 less than absolute integer midpoint distance')
    assert (len(cohort), cohort.regulated.sum(), cohort.gene.nunique()) == (3852, 101, 54)
    return cohort, excluded, qc

def ranked_ap_setup(y, scores):
    """Tie-aware AP weighted by cluster multiplicity, sorted just once."""
    order = np.argsort(-scores, kind='stable')
    ends = np.r_[np.flatnonzero(np.diff(scores[order]) != 0), len(scores)-1]
    yy = y[order].astype(float)
    def calculate(weights):
        w = weights[order]
        tp = np.cumsum(w*yy)[ends]
        all_count = np.cumsum(w)[ends]
        if not len(tp) or tp[-1] == 0:
            return np.nan
        precision = np.divide(tp, all_count, out=np.zeros_like(tp), where=all_count > 0)
        return float(np.sum(np.diff(np.r_[0, tp]) * precision) / tp[-1])
    return calculate

def compare(df, bootstrap, seed):
    y = df.regulated.to_numpy(dtype=int)
    abc, distance = df.abc_score.to_numpy(), -df.distance_bp.to_numpy(dtype=float)
    f_a, f_d = ranked_ap_setup(y, abc), ranked_ap_setup(y, distance)
    ones = np.ones(len(df))
    a, d = f_a(ones), f_d(ones)
    assert np.isclose(a, average_precision_score(y, abc), atol=1e-12)
    assert np.isclose(d, average_precision_score(y, distance), atol=1e-12)
    genes, codes = np.unique(df.gene.to_numpy(), return_inverse=True)
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(bootstrap):
        mult = rng.multinomial(len(genes), np.repeat(1/len(genes), len(genes)))
        weights = mult[codes]
        aa, dd = f_a(weights), f_d(weights)
        if np.isfinite(aa) and np.isfinite(dd):
            draws.append([aa, dd, aa-dd])
    intervals = np.percentile(draws, [2.5, 97.5], axis=0)
    metric = dict(n=len(df), positives=int(y.sum()), genes=len(genes), prevalence=float(y.mean()), missing=int((~df.matched).sum()),
                  abc_ap=a, distance_ap=d, delta=a-d, abc_low=intervals[0,0], abc_high=intervals[1,0],
                  distance_low=intervals[0,1], distance_high=intervals[1,1], delta_low=intervals[0,2], delta_high=intervals[1,2],
                  bootstrap_draws=bootstrap, bootstrap_valid=len(draws), bootstrap_seed=seed)
    curves = {}
    for name, scores in [('abc', abc), ('distance', distance)]:
        precision, recall, _ = precision_recall_curve(y, scores)
        curves[name] = [[float(r), float(p)] for r,p in zip(recall[::-1], precision[::-1])]
    return metric, curves

def verify_official(cohort, results):
    source = results / 'official/processed_merged.tsv.gz'
    if not source.exists():
        raise RuntimeError('Missing official numerical table. Run Rscript scripts/reproduce_official.R first (or retain checked-in results/official).')
    original = pd.read_csv(source, sep='\t')
    ours = cohort.set_index('id')
    for predictor, column in [('ABC.powerlaw.Score', 'abc_score'), ('baseline.distToTSS', 'distance_bp')]:
        other = original[original.pred_uid == predictor].set_index('name')
        assert set(other.index) == set(ours.index), f'Official pair set differs: {predictor}'
        other = other.loc[ours.index]
        np.testing.assert_array_equal(other.Regulated.to_numpy(), ours.regulated.to_numpy())
        np.testing.assert_allclose(other.pred_value.to_numpy(), ours[column].to_numpy(), rtol=1e-12, atol=1e-12)
        if predictor.startswith('ABC'):
            np.testing.assert_array_equal(other.Prediction.to_numpy().astype(bool), ours.matched.to_numpy())
    return {'status': 'PASS', 'pair_ids_labels_scores_distance_coverage': 'All 3,852 pairs agree with original R functions',
            'table_sha256': hashlib.sha256(source.read_bytes()).hexdigest()}

def select_examples(df):
    """Explicit illustrative rules; not a random or independent validation set."""
    definitions = [
        ('Distal agreement', df[(df.distance_bp > 100000) & df.regulated & (df.abc_score >= .02)], 'abc_score', False,
         'Highest ABC score among distal (>100 kb) positive pairs above 0.02; distance alone ranks this lower.'),
        ('Strong score, no down call', df[(~df.regulated) & (~df.significant) & (df.abc_score >= .02)], 'abc_score', False,
         'Highest ABC score among nonsignificant pairs above 0.02; a prediction–assay disagreement, not proof of no regulatory effect.'),
        ('Nearby missed effect', df[df.regulated & df.matched & (df.abc_score < .02)], 'distance_bp', True,
         'Nearest covered positive pair below the supplied 0.02 threshold; short distance can retain a useful candidate.'),
        ('Coverage gap', df[df.regulated & (~df.matched)], 'padj', True,
         'Lowest adjusted p-value among positive pairs with no matching saved prediction; zero is an evaluation fill, not a biological prediction.'),
        ('Opposite effect direction', df[(~df.regulated) & df.significant & (df.effect > 0)], 'abc_score', False,
         'Highest ABC score among significant increases; the benchmark target is down-regulation, not any expression change.')]
    examples = []
    for title, sub, sort, ascending, reason in definitions:
        assert len(sub)
        row = sub.sort_values([sort, 'id'], ascending=[ascending, True]).iloc[0]
        examples.append(dict(id=row.id, title=title, reason=reason))
    return examples

def make_figures(df, metrics, curves, examples, directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False,
                                'axes.spines.right': False, 'svg.fonttype': 'none', 'svg.hashsalt': 'engreitz-20261029'})
    directory.mkdir(exist_ok=True, parents=True)
    teal, gray = '#1f7a8c', '#555555'
    def save(fig, name):
        fig.tight_layout()
        for ext in ('svg', 'png'):
            fig.savefig(directory / f'{name}.{ext}', dpi=170, metadata={'Date': None} if ext == 'svg' else {})
        plt.close(fig)
    m = metrics[0]
    fig, ax = plt.subplots(figsize=(7.2,4.7))
    for name, label, color, style in [('abc','Saved ABC (power law)',teal,'-'), ('distance','Distance to TSS',gray,'--')]:
        c = np.array(curves['all|all'][name])
        ax.step(c[:,0], c[:,1], where='pre', color=color, ls=style, lw=1.8, label=f'{label} · AP {m[name+"_ap"]:.3f}')
    ax.axhline(m['prevalence'], color='#999', ls=':', label=f'Positive fraction {m["prevalence"]:.1%}')
    ax.set(xlim=(0,1), ylim=(0,1.03), xlabel='Recall · fraction of 101 down-regulating pairs recovered', ylabel='Precision · fraction of selected pairs down-regulating', title='Same 3,852 tested pairs; no untested negatives')
    ax.legend(frameon=False, loc='upper right', fontsize=9)
    save(fig, '01-pr')
    fig, ax = plt.subplots(figsize=(7.2,4.8))
    bins = [m for m in metrics if m['cohort']=='all' and m['range']!='all']
    for offset, name, label, color, marker in [(-.10,'abc','Saved ABC',teal,'o'),(.10,'distance','Distance',gray,'s')]:
        vals = np.array([m[name+'_ap'] for m in bins])
        lo = np.array([m[name+'_low'] for m in bins]); hi=np.array([m[name+'_high'] for m in bins])
        ax.errorbar(np.arange(3)+offset, vals, yerr=[vals-lo,hi-vals], fmt=marker, color=color, capsize=3, label=label)
    for i, item in enumerate(bins):
        ax.plot([i-.28,i+.28],[item['prevalence']]*2, color='#aaa', ls=':')
    labels=[f'{m["range"]}\n{m["n"]:,} pairs · {m["positives"]} positives\npositive fraction {m["prevalence"]:.2%}' for m in bins]
    ax.set(xticks=np.arange(3), xticklabels=labels, ylabel='Average precision (AP)', ylim=(0,1.03), title='Distance strata · gene-cluster bootstrap 95% intervals')
    ax.legend(frameon=False)
    ax.text(.01,.02,'Dotted lines: positive fraction. Compare methods within each bin.',transform=ax.transAxes,fontsize=8)
    save(fig, '02-distance')
    fig, ax = plt.subplots(figsize=(7.2,4.7))
    for positive, color, label in [(False,'#aaa','No down-regulation call'), (True,teal,'Down-regulation call')]:
        sub=df[df.regulated==positive]
        ax.scatter(sub.distance_bp/1000, sub.abc_score, s=13 if positive else 8, c=color, alpha=.75 if positive else .22, label=label)
    missing=df[~df.matched]
    ax.scatter(missing.distance_bp/1000, missing.abc_score, marker='x', c='#555', s=12, alpha=.35, label='No matched prediction (filled 0)')
    ax.axhline(.02, color='#777',ls='--',lw=1,label='Supplied ABC cutoff 0.02')
    for i, example in enumerate(examples):
        r=df[df.id==example['id']].iloc[0]
        ax.scatter([r.distance_bp/1000],[r.abc_score],s=75,facecolors='none',edgecolors='#222',lw=1)
        ax.annotate(str(i+1), (r.distance_bp/1000,r.abc_score), xytext=(5,7),textcoords='offset points',fontsize=10)
    ax.set_xscale('log'); ax.set_yscale('symlog',linthresh=.001)
    ax.set_ylim(-.00015, float(df.abc_score.max())*1.8)
    ax.set(xlabel='Distance to TSS (kb, log scale)',ylabel='Summed saved ABC score (symlog near zero)')
    ax.set_title('Agreement and disagreement · numbered cases in the explorer',pad=52)
    ax.legend(frameon=False,fontsize=8,loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2)
    save(fig,'03-disagreement')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, help='Existing checksum-verified upstream checkout or input directory')
    parser.add_argument('--bootstrap', type=int, default=2000)
    args=parser.parse_args()
    assert args.bootstrap >= 100, 'Use at least 100 resamples; release uses 2,000.'
    inputs=fetch(args.inputs)
    results=ROOT/'results'; results.mkdir(exist_ok=True)
    df, excluded, qc=build_cohort(inputs)
    qc['official_comparison']=verify_official(df, results)
    metrics=[]; curves={}
    for cohort_idx, policy in enumerate(['all','matched']):
        base=df if policy=='all' else df[df.matched]
        for bin_idx,(name,lo,hi) in enumerate(BINS):
            sub=base[(base.distance_bp>=lo)&(base.distance_bp<hi)]
            m,c=compare(sub,args.bootstrap,SEED+cohort_idx*10+bin_idx)
            m.update(cohort=policy,range=name)
            metrics.append(m); curves[policy+'|'+name]=c
    # Conventional midpoint units change at most 1.5 bp, checked as a sensitivity.
    qc['absolute_midpoint_distance_ap']=float(average_precision_score(df.regulated,-df.absolute_midpoint_distance_bp))
    qc['midpoint_distance_max_difference_bp']=float((df.absolute_midpoint_distance_bp-df.distance_bp).abs().max())
    threshold = df.abc_score >= .02
    qc['abc_threshold_counts']={'true_positive':int((threshold&df.regulated).sum()),'not_down_called':int((threshold&~df.regulated).sum()),
                                'missed_positive':int((~threshold&df.regulated).sum()),'low_score_no_down_call':int((~threshold&~df.regulated).sum())}
    examples=select_examples(df)
    df.to_csv(results/'cohort.csv',index=False, float_format='%.15g')
    excluded.to_csv(results/'excluded_rows.csv',index=False)
    pd.DataFrame(metrics).to_csv(results/'metrics.csv',index=False, float_format='%.15g')
    pd.DataFrame(examples).merge(df,on='id').to_csv(results/'examples.csv',index=False, float_format='%.15g')
    longcurves=[dict(cohort=key.split('|')[0],range=key.split('|')[1],method=method,recall=r,precision=p) for key,c in curves.items() for method,points in c.items() for r,p in points]
    pd.DataFrame(longcurves).to_csv(results/'precision_recall.csv',index=False,float_format='%.15g')
    provenance=json.loads((ROOT/'data/input_manifest.json').read_text())
    official_table=pd.read_csv(results/'official/performance_summary.tsv',sep='\t').set_index('pred_uid')
    bundle=dict(metrics=metrics,curves=curves,pairs=json.loads(df.to_json(orient='records',double_precision=15)),examples=examples,qc=qc,
                official={'abc_auprc':float(official_table.loc['ABC.powerlaw.Score','AUPRC']),
                          'distance_auprc':float(official_table.loc['baseline.distToTSS','AUPRC'])},provenance=provenance)
    (ROOT/'data.js').write_text('window.BENCHMARK = '+json.dumps(bundle,ensure_ascii=False,allow_nan=False,separators=(',',':'))+';\n')
    (results/'qc.json').write_text(json.dumps(qc,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    (results/'summary.json').write_text(json.dumps({k:v for k,v in bundle.items() if k not in ['pairs','curves']},indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    (results/'environment.json').write_text(json.dumps(dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,
                                               sklearn=__import__('sklearn').__version__,matplotlib=__import__('matplotlib').__version__,seed=SEED),indent=2)+'\n')
    make_figures(df,metrics,curves,examples,ROOT/'figures')
    print(json.dumps({'qc':qc,'metrics':metrics,'examples':examples},indent=2,ensure_ascii=False))

if __name__=='__main__':
    main()
