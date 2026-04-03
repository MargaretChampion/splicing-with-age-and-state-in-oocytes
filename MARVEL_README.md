Fixes applied to make MARVEL SE quantification work:

1. tran_id format from Preprocess_rMATS was incompatible with MARVEL's literal split on ":+@" / ":-@".
   → duplicated strand before '@'

2. SpliceJunction rownames did not match expected "chr:start:end" format
   → rebuilt coord.intron and added "chr" prefix

3. MARVEL expects coord.intron column; added explicitly

Result:
- 18,995 SE events input
- 7,552 validated and quantified
