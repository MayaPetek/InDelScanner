# Script to perform pairwise sequence alignment using mafft, between a given read and the reference
# genome found in "reference".

# mafft is a bit annoying and needs us to write an input file like this:
# > read
# LETTERSANDSHIT
# > reference
# DIFFERENTLETTERSANDSHIT

# Inlined the reference again for perf (saved file io :P)
printf "> read\n$1\n> reference\nATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCCGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTCTGACGTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCACATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCCTTCAAAGATGACGGGACCTACAAGACGCGTGCCGAGGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGTCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAAAAGCTT" > /tmp/mafftIn

# These mafft options are described by the manual as "probably most accurate, but very slow".
# `--nofft` is passed because:
#
#       --nofft
#           Do not use FFT approximation in group-to-group alignment.  Default: off
#
# "approximation" is not an encouraging thing to have turned on, and we're in no need of extra speed.
# This might yield very slightly better diffs in some edgecases, then.
# You may want to fiddle about with the various tuning options (none of which I understand). The
# defaults produce sane-looking results, but they might be useful. See `man mafft`.
#
# It's still plenty fast enough for us, since our sequence is short and we're not trying to
# match zillions of sequences together.
mafft --nofft --reorder --quiet --localpair --maxiterate 1000 /tmp/mafftIn 2>/dev/null | awk '{print toupper($0)}'
