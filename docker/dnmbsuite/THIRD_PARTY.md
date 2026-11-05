# Native tool provenance

The base image is `ghcr.io/jaeyoonsung/dnmbsuite@sha256:2702a0615bfadd6e98b8ac67b38f45d5db6cbe4957aa90a9d371e9ab62d1f47c` (Linux amd64). Its OCI source revision is `f028b1dfc1f5d15b678fa4f69dcd2775665761a3` in <https://github.com/JAEYOONSUNG/DNMBsuite>. Base tools retain their individual licenses.

`install_native_runtime.py` installs public MIT-licensed sources into a separate R library and keeps exact source and installation receipts:

- DNMB: <https://github.com/JAEYOONSUNG/DNMB/tree/976489e9967e9e4dfacc61f7e83e66b1da0f4f2b>
- Historical R DefenseViz helper namespace: <https://github.com/JAEYOONSUNG/DefenseViz/tree/53a1dc4b8c6595945762e21354fbdc86d2b3731a>

CATH `cath-resolve-hits` v0.16.10 is GPL-3.0; its license is included as `CATH-LICENSE.txt`. The Linux x86-64 static executable is downloaded from the official release and verified against SHA256 `538a43d60df049a4959cb8a2ca8d704a21b6dbd10c439384aedbdb18bd02b5b7`.

- Release: <https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.16.10>
- Exact corresponding source: <https://github.com/UCLOrengoGroup/cath-tools/tree/99edb286b0ee3a339723264486191f9f2817607d>
- Source archive: <https://github.com/UCLOrengoGroup/cath-tools/archive/99edb286b0ee3a339723264486191f9f2817607d.tar.gz>

The recipe does not redistribute REBASE/HMM databases, sample genomes or proprietary fonts. Supply database resources and any licensed Arial fonts separately. A published image must retain these notices and the upstream source/license information; the project MIT license does not replace the licenses of bundled tools.
