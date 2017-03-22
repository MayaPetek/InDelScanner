for i in $(cat alignments); do
    ./compareAcids.py $(cat reference) $i
    echo "---"
done
