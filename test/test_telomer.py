


from telofinder import get_data
from telofinder import analyze_telom_length  as alt
filename = get_data("test.telo.blocks.fasta")

def test_telomer():
    # put your code here below and remove the pass statement
    df = alt.run_telofinder(filename)
