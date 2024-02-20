# project 2: qiime2 steps

## importing

Human:

```bash
qiime tools import \
--type "SampleData[SequencesWithQuality]" \
--input-format SingleEndFastqManifestPhred33V2 \
--input-path /home/qiime2/data/metadata/human_manifest.tsv \
--output-path demux_human.qza
```

Visualizing demux file:

```bash
qiime demux summarize \
--i-data demux_human.qza \
--o-visualization demux_human.qzv
```

Mouse:

## denoising

Human:

Mouse:
