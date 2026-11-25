# Generation of Anand2013.txt and McElvain2022.txt

The nuclear responses of the models `Anand2013` and `McElvain2022` are the exact polynomials of
[DMFormFactor v6](https://www.ocf.berkeley.edu/~nanand/software/dmformfactor/) (Anand, Fitzpatrick, Haxton,
arXiv:1308.6288): `ResponseNuclear[y, j, t1, t2] = exp(-2y) sum_k a_k y^k`, with j = 1 ... 8 (M, Sigma'',
Sigma', Phi'', Phi-tilde', Delta, M-Phi'', Sigma'-Delta) and t1, t2 = 0, 1 (isospin +, -).

- `Anand2013`: the default one-body density matrices (OBDMs) of DMFormFactor,
  `SetIsotope[Z, A, "default", "default"]`.
- `McElvain2022`: the OBDMs of the [Elastic repository](https://github.com/Berkeley-Electroweak-Physics/Elastic)
  (arXiv:2208.07945), one interaction per nucleus (list in `export_responses.wls` and in the header of
  `McElvain2022.txt`). `SetIsotope[Z, A, "default", "default"]` sets b, J and T = |N - Z|/2 (for C, O, Al, Ar, Ca
  and Cs it prints "Cannot find default matrix", which is harmless), then the rows of the Elastic file replace
  the density matrix of DMFormFactor. The Elastic files (blocks with the header `2Ji 2Ti 2Jf 2Tf 2Jop 2Top` and
  rows `N 2j N' 2j' value`) use the convention of the internal tables `{value, Jop, Top, N, 2j, N', 2j'}` of
  DMFormFactor, but its `ReadInDMFile` does not read them. The script checks that J and T of the Elastic file are
  those of DMFormFactor.

## Usage

Neither DMFormFactor nor the Elastic OBDMs are copied into RAPIDD: pass their paths.

    W="/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/wolframscript"
    "$W" -file export_responses.wls <path>/dmformfactor-V6.m <path>/Elastic Anand2013 export_Anand2013.txt
    "$W" -file export_responses.wls <path>/dmformfactor-V6.m <path>/Elastic McElvain2022 export_McElvain2022.txt
    python write_responses.py export_Anand2013.txt Anand2013 ../Anand2013.txt
    python write_responses.py export_McElvain2022.txt McElvain2022 ../McElvain2022.txt \
        --keep ../McElvain2022.txt --keep-A 128,129,130,131,132,134,136

`dmformfactor-V6.m` is the original, unmodified package; `<path>/Elastic` is a clone of the Elastic repository.
An optional sixth argument of `export_responses.wls` (comma-separated A, e.g. `19,23`) exports only some nuclei.

`write_responses.py` writes the rows in the order M, Sigma'', Sigma', Phi'', Phi-tilde', Delta, M-Phi'',
Sigma'-Delta and ++, --, +-, -+, omitting the rows whose coefficients are all 0.

The Xe rows of `McElvain2022.txt` are kept as they are (`--keep`): they come from `dmformfactor-V7.wl`
(DMFormFactor with the Elastic GCN5082 OBDMs pasted in place of its `sdgXe*` tables). `write_responses.py` prints
how much they differ from the regenerated ones (at most 2e-5 of the largest coefficient of a row). Without
`--keep` the file is entirely regenerated from the Elastic files.
