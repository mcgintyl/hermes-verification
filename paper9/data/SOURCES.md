# Data sources

The files in `famaey/` and `mistele/` are unmodified copies of files published on Zenodo under the Creative Commons Attribution 4.0 licence (CC BY 4.0). Their MD5 checksums were matched against the Zenodo record listings on 2026-09-11. `rmax_x.csv`, `ages_133.csv` and `sparc_manifest.csv` were produced by this work and are described below. `sparc/` is empty until `fetch_sparc.py` downloads the rotation curves; those are not redistributed here.

## Famaey, Pizzuti & Saltas

Famaey, B., Pizzuti, L. & Saltas, I. (2025). *On the nature of the missing mass of galaxy clusters in MOND: the view from gravitational lensing* (data and code). Zenodo. https://doi.org/10.5281/zenodo.15299349. Paper: arXiv:2410.02612.

- `famaey/ClusterInfo.py`: the gas-density fit parameters (mostly from Laudato, Salzano & Umetsu 2022, MNRAS 511, 1878, fitted to the Chandra data of Donahue et al. 2014), BCG masses, r₂₀₀, redshifts and dynamical-state flags, and the gas-mass integral.
- `famaey/fgasMACS.txt`: the member-galaxy template.
- `famaey/profiles/results_profile_NOgnfw_all_*.txt`: nine files from the record's `mass_profiles.zip` (MD5 `931bf31b3d68a04f6f436bd941549ab3`). Only the radius column is used, as the 50-point radial grid of the baryonic carrier.

## Mistele et al.

Mistele, T., Lelli, F., McGaugh, S., Schombert, J. & Famaey, B. (2025). *Mass models of galaxy clusters from a non-parametric weak-lensing reconstruction* (data). Zenodo. https://doi.org/10.5281/zenodo.15476959. Paper: arXiv:2506.13716.

- `mistele/<cluster>.csv`: the deprojected, three-dimensional enclosed mass M(<r) with its uncertainties. The code uses the `M_Msun` and `M_stat_err_Msun` columns.
- `mistele/<cluster>-M-corr.csv`: the node-to-node correlation matrix for M(<r).

## SPARC rotation curves (fetched, not redistributed)

Lelli, F., McGaugh, S. S. & Schombert, J. M. (2016). *SPARC: Mass models for 175 disk galaxies with Spitzer photometry and accurate rotation curves*. AJ 152, 157. Data: http://astroweb.case.edu/SPARC/.

`fetch_sparc.py` downloads the public archive `Rotmod_LTG.zip` (SHA-256 `0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588`, 110,737 bytes), extracts only the 133 `*_rotmod.dat` files the paper scores into `sparc/`, and verifies each one against `sparc_manifest.csv`. SPARC's terms ask that users of the data cite the paper above. If the archive is ever repackaged and its own checksum changes, the per-file checks are what matter: they confirm the rotation curves are the same bytes the published result used.

## Files produced by this work

- `rmax_x.csv`: the radial extent of the X-ray data, R_max,X, transcribed from Mistele et al. (arXiv:2506.13716v2), Table 1. It is used only to count the nodes that sit on extrapolated gas, and for the 1/r⁴ alternative gas treatment.
- `ages_133.csv`: the median stellar age t₅₀ and the 98th-percentile baryonic acceleration g₉₈ for each of the 133 galaxies, as used in Paper 1 and carried unchanged into Paper 9.
- `sparc_manifest.csv`: for each of the 133 galaxies, the rotation-curve file name, its SHA-256 and its size in bytes, taken from the official archive above.

| file | SHA-256 | bytes |
|---|---|---|
| rmax_x.csv | b2742a29a8a37cd7d2025ae3110ad2134af562f6cd2766c22a7c25c0896c0d52 | 357 |
| ages_133.csv | 9ad89b62e29f64b08c30117669c4f4a927c68e547a01cf98cc1cadb342744fd9 | 4349 |
| sparc_manifest.csv | 88f5e134f796524243839e78bba544d3f4f3e2b207af6de659d0f9e7f09c4c84 | 13075 |

The reference values the two scripts check against, in `expected/`:

| file | SHA-256 | bytes |
|---|---|---|
| expected/table2_seven.csv | 47ddd6b7867018f8b42b55a422d768205381384471c0240f918eb285e16bb96c | 647 |
| expected/table2_full_precision.csv | 76e5b9e4e6a5b5f1d6861543a273f62743c3aebac062ccc3df5483fe781f127a | 727 |
| expected/paper1_per_galaxy.csv | 622a9db40b7309b03ab50ee6e48deba42b97aedc0d60acb3cf56cbd1679315cf | 7097 |

## Checksums

| file | SHA-256 | MD5 | bytes |
|---|---|---|---|
| famaey/ClusterInfo.py | 59bfc1d0d4dc35bf9dcc3e6b02c42d4b57bbba081d6749e0b1f08b04f23f9db4 | d388b00fd22c0a2c604de792d99c9ee0 | 6969 |
| famaey/fgasMACS.txt | 7ce1d9c75cd7650c1ee8170d4ef732eba30bb2f4fecb074da66309b279a297c8 | 569b57087963c619cb39f82154bf056b | 5100 |
| famaey/profiles/results_profile_NOgnfw_all_a2261.txt | ed5513ae1d40864986524e8dadebe503c9252b701d133dadc554d373f93db1b9 | 914a2882c7aef47ee216eaabbc1006c0 | 6133 |
| famaey/profiles/results_profile_NOgnfw_all_a611.txt | e29fb674a711d10ef84ffa406f1434bdaa03955876bcbf014623b07ee1ec44ba | 4a9ba0e9f72de7d9ac0c5aab55b924d6 | 6176 |
| famaey/profiles/results_profile_NOgnfw_all_macsj0429.txt | f3d13362a5cca50e57b7c4222102af686dee9d4fda4b3932f7754be084b05e75 | a057eb6f7ec19693e8766a883b5897e2 | 6179 |
| famaey/profiles/results_profile_NOgnfw_all_macsj0744.txt | d4ab9865bdb79281393be92e48f4e0cdb7ccef88423f3b5649eb72c4ab202457 | c166395e0ea9b545e76f5c0ea492c9d1 | 6154 |
| famaey/profiles/results_profile_NOgnfw_all_macsj1115.txt | 5a0ac395e452a58175bdb2c57a855770ebef37123852e9e24983e67699666c26 | 69bec100ea620c4cf05d737b099c75c2 | 6150 |
| famaey/profiles/results_profile_NOgnfw_all_macsj1206.txt | c5a821891af2bd70ec84e003236501b9759fd73d67de3166e02403b8b9d6d0c8 | 3863dc15a6a6c19ee5912fa89514b521 | 6140 |
| famaey/profiles/results_profile_NOgnfw_all_macsj1720.txt | b3d5a046e08f792a0dde34e4e3beb743a6628927ba160d0110bcca9c2df7a6df | 53af29722ae91cb669623a73d26dfce7 | 6186 |
| famaey/profiles/results_profile_NOgnfw_all_rxj1347.txt | 37fb300640a7844b052afbe7f5c502306444934c1421ddda679b3dda60f7b486 | 328af24d9293c48a06235a8442a58980 | 6183 |
| famaey/profiles/results_profile_NOgnfw_all_rxj2129.txt | 885db0f3a5bdd29664622e1f2991291d4783aade782b6373ceb7de0b8c83a8d3 | f2702919c30083ec6488d2e974872a40 | 6170 |
| mistele/abell2261-M-corr.csv | 3ee59499072f0a7948483862f84e0fd93d82aaf7ec2c131dc98791b49261f271 | 163e44f0a1f2b9503a75b708d77e38ec | 1932 |
| mistele/abell2261.csv | 38b23772ffcec5f0cc17d925e264cb57c825a39ef2f0bf7c26f7841290cc91e1 | f01c441ae4513756af61a5c1cd730d4b | 2990 |
| mistele/abell611-M-corr.csv | 7372a8ec6d1da1298b355afadfd546de1decd949841ef84a230f9f97d69bd0f6 | bbb83af0251e018ccfdfe1507c89000f | 1963 |
| mistele/abell611.csv | ccadf240b730ca5b8ff784d292ddaf3bca2ad45df267c33dc634f22bca38310c | df8f0ce5801243396c77c1a0d8dea253 | 2991 |
| mistele/macsj0429-M-corr.csv | 8bff258c94977108854fed8691355bf9d9e58a4626bab96446cb0f174527e116 | eec89e04b75f72940f244b0270a06f17 | 1897 |
| mistele/macsj0429.csv | 97062cffe4fd31f4c94f45569471b1bd5ef349cd3dfefb7881468b6a0b8d8bf2 | 6e833cebfbc5ca71dd1a8337a4f1f362 | 2989 |
| mistele/macsj0744-M-corr.csv | 2d9b899a6b77dea29c043cbdf5be498bb9e16b601cde4d1ffc2b8f1f75f35691 | a6fd79c2bd74cf15ca3dfcab971c79e2 | 1926 |
| mistele/macsj0744.csv | 1ecb8c079209015ae36b5f3729425458a82216d294cb4ce5e8f1b1b98e96da5d | 338f99f376e1dac4fb33f47a2465fd6c | 2992 |
| mistele/macsj1115-M-corr.csv | d75d9c643f8f08724ddaff5b40d3498025d56f7b0af0ea202cbebfc2ad6f7cee | 5c17fc4d89a386f2f9e76cb1c1560ae0 | 1939 |
| mistele/macsj1115.csv | 301b9d75b7e9c11bacf736d2a7f5a4f05a62ef2eb798df672c2433bb98cf4cce | 7772145ca1b2d989ec489de96982b8f8 | 2998 |
| mistele/macsj1206-M-corr.csv | b391770ade51da319061408f9eb0e096ea511c13edcb086d1b757685e2138bde | 33a30dc40b74de73310416699042feff | 1975 |
| mistele/macsj1206.csv | 84f72e560030556429f39e96e08c185e2c8352e9c4a98634931c2beb8d9158a1 | 406da1b0be1230fee4f7d87cfb588e40 | 3005 |
| mistele/macsj1720-M-corr.csv | ad1ae0bbdfb5ce1d7e4108e6eac29ad40f9d0b3808e070325ff8b6faa00415b5 | ff28190c981ce04fd2a1f7457e2e43f0 | 1929 |
| mistele/macsj1720.csv | 4bf85ced0f147ba91cfbeac5d108cf1df83f3d5d0863ae6ff43b6e570cd73edb | 27c36e096e3e0bd74111f172bb7b8523 | 2990 |
| mistele/rxj1347-M-corr.csv | ca6491730ce09c99789509491471161547e866ebc30d3b826bae46834b82a851 | 8c9bfd485596743fe25a4664af43ffb5 | 1934 |
| mistele/rxj1347.csv | 1e3e5dc6a3ed63946c4e934715eb374aafe30cf7390b08423ddfb1e23d5c83ee | c656e4d7ae01d64392906ea16bc860e8 | 2995 |
| mistele/rxj2129-M-corr.csv | 81b04f67cf5e84deba6075eae779260dcaa2d8e000792ae7838004196a95755f | 3f4ca3f56d76a1d9c330e694cf11f17c | 1968 |
| mistele/rxj2129.csv | 4828b3c83bcab6ada0314cc680e265d234d4b0eb1ac3dae514ea5e4f029b8908 | ece2b47aee9f14f15b3e7c3a718ee48f | 2996 |

The SHA-256 of every one of the 133 SPARC rotation curves is in `sparc_manifest.csv`.
