"""Journal display export of audited figures; source numbers remain byte-exact."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone
import matplotlib.pyplot as plt
import pandas as pd
import plot_rpe1_measurement as original

BASE = Path(__file__).resolve().parent / 'results/rpe1_measurement_validation_v1'
OUT = BASE / 'journal_figures'
TITLES = {
    'Figure_1_baseline': 'Baseline response utility and target identity',
    'Figure_2_factorial_q025': 'Separating model input and covariate measurement effects',
    'Figure_3_linear_identity': 'Linear model identity across measurement conditions',
    'Figure_S1_covariate_criteria': 'Covariate measurement contrasts and operational references',
    'Figure_S2_signed_prediction_MSE': 'Signed prediction error across model input measurements',
}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def journal_save(fig, out, stem, caption):
    # Overall titles and explanatory captions belong to manuscript legends.
    for text in list(fig.texts):
        text.remove()
    if stem == 'Figure_S1_covariate_criteria':
        caption += (' Dashed references are mean agreement ≥0.90, target p10 ≥0.80, and scaled absolute utility change ≤0.10. '
                    'Cross symbols and asterisks mark native all-floor cases; daggers mark at least one thinning seed crossing the ranking criterion. '
                    'The training-mean/Dc64/fixed Sq mean is at the numerical 0.90 boundary; its actual seed crossings remain recorded.')
    sizes = {
        'Figure_1_baseline': (170, 122), 'Figure_2_factorial_q025': (170, 180),
        'Figure_3_linear_identity': (170, 158), 'Figure_S1_covariate_criteria': (170, 168),
        'Figure_S2_signed_prediction_MSE': (170, 100),
    }
    width, height = sizes[stem]
    fig.set_size_inches(width / 25.4, height / 25.4)
    for ax in fig.axes:
        for txt in [ax.title, ax._left_title, ax._right_title, ax.xaxis.label, ax.yaxis.label, *ax.texts, *ax.get_xticklabels(), *ax.get_yticklabels()]:
            txt.set_fontsize(7.5)
        ax.tick_params(labelsize=7)
    for legend in fig.legends:
        for text in legend.get_texts():
            text.set_fontsize(7.5)
            text.set_text(text.get_text().replace('Small open symbols: fold means', 'Open symbols: folds').replace('Five permutation draws and range', 'Five permutations and range'))
    if stem == 'Figure_1_baseline':
        fig.subplots_adjust(left=.29, right=.985, top=.80, bottom=.14, wspace=.45)
        fig.axes[0].set_xlabel('Total utility\n(×10⁻³)', fontsize=7.5)
        # The empirical-anchor fold reaches 101.295; retain the entire point.
        fig.axes[0].set_xlim(-3, 110)
        fig.axes[1].set_xlabel('Identity gain\n(×10⁻³)', fontsize=7.5)
        fig.axes[1].set_title('B  Model identity', loc='left', fontsize=7.5)
        fig.axes[2].set_title('C  Empirical A', loc='left', fontsize=7.5)
        fig.axes[2].set_ylabel('Identity gain (×10⁻³)', fontsize=7.5)
        for text in list(fig.axes[1].texts):
            if 'Technical anchor' in text.get_text(): text.remove()
    elif stem == 'Figure_2_factorial_q025':
        fig.subplots_adjust(left=.13, right=.985, top=.88, bottom=.07, wspace=.42, hspace=.20)
        for ax, label in zip(fig.axes[:4], ['scTenifoldKnk*', 'Pearson', 'Held-out linear', 'Training mean']):
            ax.set_title(label, loc='left', fontsize=7.5)
    elif stem == 'Figure_3_linear_identity':
        fig.subplots_adjust(left=.12, right=.985, top=.86, bottom=.08, wspace=.62, hspace=.37)
        for i, ax in enumerate(fig.axes):
            ax.set_ylabel('Identity gain (×10⁻³)' if i % 3 == 0 else 'EDR − raw identity\n(×10⁻³)', fontsize=7.5)
    elif stem == 'Figure_S1_covariate_criteria':
        fig.subplots_adjust(left=.39, right=.985, top=.92, bottom=.08, wspace=.36)
        for ax, label in zip(fig.axes, ['A  Mean overlap', 'B  Target p10', 'C  Scaled change']):
            ax.set_title(label, loc='left', fontsize=7.5)
            ax.set_xlabel('Readout', fontsize=7.5)
        labels=[t.get_text().replace('Pearson association','Pearson').replace('Held-out linear','Linear').replace('scTenifoldKnk','Network') for t in fig.axes[0].get_yticklabels()]
        fig.axes[0].set_yticklabels(labels, fontsize=7)
    else:
        fig.subplots_adjust(left=.13, right=.985, top=.80, bottom=.19, wspace=.36)
        for ax, label in zip(fig.axes, ['A  Response arm B', 'B  Response arm C']):
            ax.set_title(label, loc='left', fontsize=7.5)
            ax.set_xlabel('Molecule retention q', fontsize=7.5)
    fig.savefig(out / (stem + '.pdf'))
    fig.savefig(out / (stem + '.png'), dpi=300)
    plt.close(fig)
    original.CAPTIONS[stem] = caption


def main():
    if OUT.exists(): raise RuntimeError('Refusing to overwrite journal figures')
    evaluation = BASE / 'evaluation'
    manifest = json.loads((evaluation / 'manifest.json').read_text())
    for name in original.INPUTS:
        assert sha(evaluation / name) == manifest['output_sha256'][name]
    data = {name: pd.read_csv(evaluation / name) for name in original.INPUTS}
    OUT.mkdir()
    original.save = journal_save
    original.style()
    original.fig1(data, OUT); original.fig2(data, OUT); original.fig3(data, OUT)
    original.fig_s1(data, OUT); original.fig_s2(data, OUT)
    (OUT / 'Figure_legends.md').write_text('# Journal figure titles and legends\n\n' + '\n\n'.join(
        f'## {stem}\n\n**{TITLES[stem]}.** {caption}' for stem, caption in original.CAPTIONS.items()) + '\n')
    exact = []
    for p in sorted(OUT.glob('*.csv')):
        assert sha(p) == sha(BASE / 'figures' / p.name), p.name
        exact.append(p.name)
    doc = {'status': 'journal_display_export_pending_visual_review', 'created_utc': datetime.now(timezone.utc).isoformat(),
           'source_script_sha256': sha(Path(original.__file__)), 'export_script_sha256': sha(Path(__file__)),
           'evaluation_manifest_sha256': sha(evaluation / 'manifest.json'),
           'source_csvs_exactly_unchanged': exact, 'display_only_changes': ['170 mm figure width', '300 dpi PNG and embedded vector PDF fonts', '7 to 7.5 pt labels at final width', 'Overall titles and captions in manuscript legend file', 'Shortened keys and panel labels'],
           'output_sha256': {p.name: sha(p) for p in sorted(OUT.iterdir()) if p.is_file()}}
    (OUT / 'manifest.json').write_text(json.dumps(doc, indent=2) + '\n')
    print(json.dumps({'out': str(OUT), 'exact_source_tables': len(exact)}))


if __name__ == '__main__': main()
