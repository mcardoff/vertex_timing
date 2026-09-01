#!/usr/bin/env python3
"""Splice the per-step input lists and importance charts from
results/model_inputs.json into results/method_explainer.html.

Idempotent: each block is fenced by an HTML comment marker, so re-running
replaces rather than duplicates.
"""
import json, re, sys

R = json.load(open("results/model_inputs.json"))
H = open("results/method_explainer.html").read()

TAGDER = set(R["stage1"]["tagger_derived"]) | {"ph"} | set(R["tagger_r2"]["new_features"])


def bars(rows, key="imp", unit="", tag=TAGDER, n=10):
    rows = rows[:n]
    mx = max(r[key] for r in rows) or 1.0
    out = ['<div class="imp">']
    for r in rows:
        w = 100.0 * r[key] / mx
        cls = ' class="t"' if r["f"] in tag else ""
        out.append(f'<span class="f">{r["f"]}</span>'
                   f'<span class="bar"><i{cls} style="width:{w:.1f}%"></i></span>'
                   f'<span class="v">{r[key]:.1f}{unit}</span>')
    out.append("</div>")
    return "".join(out)


LEG = ('<p class="lg"><span><s class="a"></s>reconstructed input</span>'
       '<span><s class="t"></s>from the tagger &#9312;</span></p>')


def spec(*pairs):
    return ('<div class="spec">'
            + "".join(f"<span>{k} <b>{v}</b></span>" for k, v in pairs)
            + "</div>")


def flist(title, feats):
    return (f'<details><summary>{title} &mdash; {len(feats)} columns</summary>'
            f'<p class="fl">{", ".join(feats)}</p></details>')


B = {}

# ---- 1 tagger ------------------------------------------------------------
r1, r2 = R["tagger_r1"], R["tagger_r2"]
B["s1"] = (
    spec(("model", "XGBoost &#183; depth 7 &#183; 600 trees"),
         ("inputs", f"{r1['n_features']} &rarr; {r2['n_features']}"),
         ("rows", f"{r1['n_rows']:,} timed tracks"),
         ("label", "truth_is_hs"),
         ("positives", f"{r1['hs_frac']}%"))
    + "<p><strong>Round 1 &mdash; what it leans on.</strong> Gain-normalised, "
      "as a percentage of the model's total:</p>"
    + bars(r1["top"], tag=set())
    + '<div class="box w"><span class="k">One feature carries two-thirds of it</span>'
      f'<p><span class="mono">closer_to_pu_than_pv</span> alone is <strong>{r1["top"][0]["imp"]}%</strong> '
      'of round 1&rsquo;s gain. It is a <em>z-geometry</em> question &mdash; is this track nearer '
      'some pileup vertex than the primary &mdash; and it contains no timing information at all. '
      'The round-1 tagger is therefore mostly a restatement of the track-to-vertex association, '
      'which is exactly why round 2 was needed.</p></div>'
    + "<p><strong>Round 2 &mdash; after adding neighbour context:</strong></p>"
    + bars(r2["top"]) + LEG
    + f'<p>The six neighbour columns take <strong>{r2["neighbor_share"]}%</strong> of the total '
      f'between them, and <span class="mono">nb_phzt</span> &mdash; hard-scatter probability mass '
      'within 2&nbsp;mm in z <em>and</em> 60&nbsp;ps in time &mdash; is second overall at '
      f'{r2["top"][1]["imp"]}%. Note what happens to the z-geometry feature: it falls from '
      f'{r1["top"][0]["imp"]}% to {r2["top"][0]["imp"]}%. The neighbourhood evidence '
      '<em>displaces</em> the single-feature reliance rather than piling on top of it, which is '
      'the sign that round 2 added information rather than redundancy.</p>'
    + flist("Full round-1 input list", r1["features"])
    + flist("Round-2 additions", r2["new_features"]))

# ---- 2 learned score -----------------------------------------------------
s1 = R["stage1"]
B["s2"] = (
    spec(("model", "XGBoost &#183; depth 6 &#183; 500 trees &#183; 5 seeds"),
         ("inputs", str(s1["n_features"])),
         ("rows", f"{s1['n_rows']:,} candidate clusters"),
         ("label", "|t_tagw &#8722; t_truth| &lt; 60 ps"))
    + bars(s1["top"]) + LEG
    + '<div class="box a"><span class="k">The learned score is not a better &Sigma;p<sub>T</sub></span>'
      f'<p>The eleven tagger-derived columns take <strong>{s1["tagger_share"]}%</strong> of the '
      f'total gain, and the raw <span class="mono">trkptz_score</span> takes '
      f'<strong>{s1["trkptz_share"]}%</strong> &mdash; the old formula is, as a value, ignored. '
      '(Its within-event <em>normalised</em> form, '
      '<span class="mono">trkptz_ratio_to_max</span>, does survive at 6.0%: the ranking is worth '
      'something, the magnitude is not.)</p>'
      '<p>The top two features are both <em>agreement</em> measures &mdash; '
      '<span class="mono">dt_tagw_agg</span> and <span class="mono">dt_to_agg</span>, the distance '
      'from this cluster&rsquo;s time to the whole-event aggregate of &#9315;, together '
      f'{s1["top"][0]["imp"] + s1["top"][1]["imp"]:.1f}%. So the selector is mostly asking '
      '<strong>&ldquo;does this cluster&rsquo;s time agree with what all the tracks in the event say?&rdquo;</strong> '
      'That is the reason &#9315; cannot be treated as an optional extra bolted on at the end: '
      'the aggregate answer is also the reference the selector measures candidates against.</p></div>'
    + flist("Full input list", s1["features"]))

# ---- 3 re-timing ---------------------------------------------------------
a = R["step3_ablation"]
d1 = a["tagw_no_winsor"] - a["cluster_time_invvar"]
d2 = a["tagw_double_winsor"] - a["tagw_no_winsor"]
B["s3"] = (
    spec(("model", "none &mdash; closed form"),
         ("inputs", "p&#7522; and &#963;&#7522; per track"),
         ("changes", "the reported time only"))
    + '<p>Because there is no model here, the impact decomposes exactly. Holding the selection '
      'fixed at <em>TRKPTZ&rsquo;s own choice</em> of cluster and changing only the arithmetic, '
      'on the 94k-event Z+jets training set:</p>'
    + '<div class="tw"><table><thead><tr><th>time estimator on the TRKPTZ-picked cluster</th>'
      '<th class="n">core frac</th><th class="n">&#916;</th></tr></thead><tbody>'
      f'<tr><td>inverse-variance mean <span class="mono">(today)</span></td>'
      f'<td class="n">{a["cluster_time_invvar"]:.2f}</td><td class="n">&mdash;</td></tr>'
      f'<tr><td>+ p&#7522; weighting</td><td class="n">{a["tagw_no_winsor"]:.2f}</td>'
      f'<td class="n up">+{d1:.2f}</td></tr>'
      f'<tr class="hl"><td>+ double Winsorisation</td>'
      f'<td class="n"><strong>{a["tagw_double_winsor"]:.2f}</strong></td>'
      f'<td class="n up">+{d2:.2f}</td></tr>'
      '</tbody></table></div>'
    + f'<p>Nearly all of it &mdash; <strong>{d1:.2f} of the {d1 + d2:.2f} points</strong> &mdash; is the '
      'probability weighting; the two Winsorisation passes add only '
      f'{d2:.2f}. The trimming is cheap insurance against a single wild track, not the mechanism.</p>'
    + '<p class="note" style="border:0;padding:0;margin:.4rem 0 0">Measured over all 93,977 '
      'novbs-augmented events, so the 62.12 baseline here is not the 61.87 headline figure, which '
      'is the canonical-selection subset. The deltas are what transfer.</p>')

# ---- 4 e2e ---------------------------------------------------------------
e = R["e2e"]
B["s4"] = (
    spec(("model", "MLP 42&#8594;96&#8594;96&#8594;1 + 5 unrolled IRLS"),
         ("inputs", str(e["n_features"])),
         ("training", "full-batch Adam &#183; 200 epochs"),
         ("loss", "smoothed |t&#8320;&#8722;t_truth| &lt; 60 ps"))
    + '<p>Permutation importance: shuffle one input column across the held-out fold, re-run the '
      'whole weight-and-IRLS forward pass, and record how many <em>points of core fraction</em> '
      f'are lost. The unshuffled fold scores {e["fold0_core_frac"]}%.</p>'
    + bars(e["top"], key="drop", unit="") + LEG
    + '<div class="box a"><span class="k">This network is a refinement on the tagger, not a rival to it</span>'
      f'<p>Shuffling <span class="mono">ph</span> alone costs <strong>{e["top"][0]["drop"]:.1f} points</strong> '
      f'&mdash; more than a fifth of the entire answer &mdash; against {e["top"][1]["drop"]:.1f} for the next '
      f'feature and under {e["top"][3]["drop"]:.1f} for everything from fourth place down. The weights it '
      'learns are essentially p&#7522; plus a correction, with '
      '<span class="mono">timeRes</span> as the main independent contributor.</p>'
      '<p>That is the single most consequential fact on this page for anyone deciding what to '
      'reuse: <strong>the per-track tagger is the load-bearing component of all five changes.</strong> '
      'It drives &#9313;, it is two-thirds of &#9314;&rsquo;s inputs, and it is four-fifths of &#9315;.</p></div>'
    + flist("Full input list", e["features"]))

# ---- 5 meta --------------------------------------------------------------
m = R["meta"]
B["s5"] = (
    spec(("model", "XGBoost &#183; depth 4 &#183; 500 trees &#183; 3 seeds"),
         ("inputs", str(m["n_features"])),
         ("rows", f"{m['n_decisive']:,} of {m['n_events']:,} events"),
         ("label", "aggregate right &amp; cluster wrong"))
    + f'<p>The chooser only ever sees the events where the two answers disagree &mdash; '
      f'<strong>{m["decisive_frac"]}%</strong> of the sample. On the other {100 - m["decisive_frac"]:.1f}% '
      'the two agree and the switch is irrelevant. Within the decisive set the aggregate answer is '
      f'the right one {m["t0_wins_frac"]}% of the time, so the default must be to keep the cluster '
      'answer and switch only on evidence.</p>'
    + bars(m["top"], tag=set())
    + '<p>Every feature at the top is a measure of <em>agreement or corroboration</em>: how far the '
      'two answers sit apart (<span class="mono">abs_dt_base_t0</span>), how confident the selector '
      'was and by how much it beat its runner-up (<span class="mono">p_gap</span>), and how many '
      'clusters and tracks independently sit near each candidate '
      '(<span class="mono">maxp_near_t0</span>, <span class="mono">n_cl_near_t0</span>). No single '
      'feature exceeds 8.2% &mdash; unlike every other model here, this one spreads its decision '
      'thinly, which is consistent with it being the weakest link in the chain.</p>'
    + flist("Full input list", m["features"]))

# ---- splice --------------------------------------------------------------
ANCH = {
    "s1": ("<span class=\"worth\">worth +1.7 pts on the Z+jets t&#8320; alone</span></p>", "after"),
    "s2": ("<span class=\"worth\">+3.0 to +5.9 pts across samples</span></p>", "after"),
    "s3": ("<h3><span class=\"num\">3</span> Re-timing the chosen cluster</h3>", "keep"),
    "s4": ("<h3><span class=\"num\">4</span> A second answer, computed from all tracks</h3>", "keep"),
    "s5": ("<h3><span class=\"num\">5</span> A chooser between the two answers</h3>", "keep"),
}

for key in ("s1", "s2", "s3", "s4", "s5"):
    blk = f"<!--IMP:{key}-->{B[key]}<!--/IMP:{key}-->"
    fence = re.compile(f"<!--IMP:{key}-->.*?<!--/IMP:{key}-->", re.S)
    if fence.search(H):
        H = fence.sub(blk, H)
        continue
    if key == "s3":
        # after the "needs no new selection" box, which closes the section
        a2 = "it can be dropped into the existing TRKPTZ pipeline on its own.</p>\n  </div>"
        assert a2 in H, key
        H = H.replace(a2, a2 + "\n  " + blk, 1)
    elif key == "s4":
        a2 = "why aggregating across the whole event can.</p>\n  </div>"
        assert a2 in H, key
        H = H.replace(a2, a2 + "\n  " + blk, 1)
    elif key == "s5":
        a2 = ("<p>Its inputs include corroboration features: how many tracks and clusters sit "
              "near each candidate answer, and how far apart the two answers are.</p>")
        assert a2 in H, key
        H = H.replace(a2, blk, 1)
    else:
        a2 = ANCH[key][0]
        assert a2 in H, key
        H = H.replace(a2, a2 + "\n  " + blk, 1)

open("results/method_explainer.html", "w").write(H)
print("spliced 5 blocks")
