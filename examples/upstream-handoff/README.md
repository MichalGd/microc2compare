# Upstream sample-sheet example

[samplesheet.csv](samplesheet.csv) uses the established microc2tracks column names. It illustrates three independent control cultures and three independent treated cultures. C1 has two independently prepared libraries; those two technical sources become one consolidated downstream observation. The sample labels alone do not establish the biological/library relationships described in this example.

FASTQ paths are placeholders and no real or synthetic contact data accompany this sheet. For your own run use its actual sheet and results. The [handoff guide](../../docs/13-upstream-handoff.md) shows selection, commands, expected draft outputs, missing-evidence handling and the different requirements for lanes of one library. Runnable numerical data are generated separately by the [tiny demonstration](../tiny/README.md).

Tests build small original Cooler files in the real upstream directory/manifest layout. They verify selection, missing/invalid evidence, draft export and unchanged upstream hashes; they do not establish production biological validity.
