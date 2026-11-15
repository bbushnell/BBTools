package hiseq;

import parse.LineParser;
import parse.LineParserS3;
import structures.ByteBuilder;

/**
 * Indexed parsing of seven-field Illumina sequencing read headers. Each
 * mutable instance retains its input and delimiter scratch for one caller.
 * Accessors parse on demand; parse() itself does not validate numeric fields.
 * Comment accessors require the corresponding space/tab-separated terms.
 *
 * Handles header formats like:
 * - @VP2-06:112:H7LNDMCVY:2:2437:14181:20134 (NovaSeq6k)
 * - @MISEQ08:172:000000000-ABYD0:1:1101:18147:1925 1:N:0:TGGATATGCGCCAATT
 * - @A00178:38:H5NYYDSXX:2:1101:3007:1000 1:N:0:CAACCTA+CTAGGTT (NovaSeq)
 * - @LH00223:28:22GLGMLT3:1:1101:5928:1016 1:N:0:CTGCTTGGTT+CTAACGACAG
 *
 * @author Brian Bushnell
 * @date April 3, 2024
 */
public class IlluminaHeaderParser2 extends ReadHeaderParser{

	/*--------------------------------------------------------------*/
	/*----------------        Expected Format       ----------------*/
	/*--------------------------------------------------------------*/

	//@VP2-06:112:H7LNDMCVY:2:2437:14181:20134 (Novaseq6k)
	//2402:6:1101:6337:2237/1
	//MISEQ08:172:000000000-ABYD0:1:1101:18147:1925 1:N:0:TGGATATGCGCCAATT
	//HISEQ07:419:HBFNEADXX:1:1101:1238:2072
	//A00178:38:H5NYYDSXX:2:1101:3007:1000 1:N:0:CAACCTA+CTAGGTT
	//@LH00223:28:22GLGMLT3:1:1101:5928:1016 1:N:0:CTGCTTGGTT+CTAACGACAG (NovaseqX)

	//	@HWI-Mxxxx or @Mxxxx - MiSeq
	//	@HWUSI - GAIIx
	//	@HWI-Dxxxx - HiSeq 2000/2500
	//	@Kxxxx - HiSeq 3000(?)/4000
	//	@Nxxxx - NextSeq 500/550
	//	@Axxxxx - NovaSeq
	//	@Vxxxxx = NextSeq 2000
	//	@AAxxxxx - NextSeq 2000 P1/P2/P3
	//	@Hxxxxxx - NovaSeq S1/S2/S4
	//
	//	AAXX = Genome Analyzer
	//	BCXX = HiSeq v1.5
	//	ACXX = HiSeq High-Output v3
	//	ANXX = HiSeq High-Output v4
	//	ADXX = HiSeq RR v1
	//	AMXX, BCXX =HiSeq RR v2
	//	ALXX = HiSeqX
	//	BGXX, AGXX = High-Output NextSeq
	//	AFXX = Mid-Output NextSeq
	//	5 letter/number = MiSeq

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();
		ihp.test(args.length>0 ? args[0] : null);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains id_ and replaces delimiter bounds; returns this mutable parser. */
	public IlluminaHeaderParser2 parse(String id_){
		id=id_;
		lp.set(id_);
		//LineParserS3 splits on ':' AND treats whitespace as a term boundary, so the comment fields
		//(pairCode/chastity/controlBits/barcode) land in their OWN terms after the coordinate block.
		//whitespaceIndex = index of the term containing/preceding the space (6 for the standard 7-coord
		//header), so pairCode is at whitespaceIndex+1, etc. This boundary behavior is what makes every
		//comment accessor below correct; without it pairCode would read the chastity term.
		whitespaceIndex=lp.indexOfWhitespace();
		return this;
	}

	/** Whether the field layout looks shrinkable and does not already use the ::: prefix. */
	public boolean canShrink(){
		return looksValid() && !looksShrunk();
	}

	/** Cheap comment-layout heuristic for shrinking, not full header or numeric validation. */
	public boolean looksValid(){
		//Shrinking requires a comment after seven coordinate fields, optionally with a third index.
		return(lp.terms()>=8 && whitespaceIndex>=6 && whitespaceIndex<=7);
	}

	/** Recognizes the leading three empty fields emitted by FastqReadInputStream's shrinker. */
	public boolean looksShrunk(){
		//Bounds are absolute exclusive ends: bounds[2]==2 identifies the leading ::: prefix.
		return(lp.terms()>3 && lp.bounds().get(2)==2);
	}

	@Override
	public String machine(){
		return lp.terms()<=0 ? null : lp.parseString(0);
	}

	@Override
	public String sample(){
		return null;
	}

	@Override
	public int run(){
		return lp.parseInt(1);
	}

	@Override
	public String flowcell(){
		return lp.parseString(2);
	}

	@Override
	public int lane(){return lp.parseInt(3);}

	/** Returns the tile number within the lane */
	@Override
	public int tile(){return lp.parseInt(4);}

	/** Returns the x-coordinate position on the tile */
	@Override
	public int xPos(){return lp.parseInt(5);}

	/** Returns the y-coordinate position on the tile */
	@Override
	public int yPos(){return lp.parseInt(6);}

	/**
	 * Returns the pair number character (1 or 2 for paired-end reads).
	 * Extracts from position after whitespace separator.
	 * @return Character indicating read pair number
	 */
	@Override
	public char pairCode(){return lp.parseChar(whitespaceIndex+1, 0);}

	/**
	 * Returns the chastity filter code (Y for failed, N for passed).
	 * Indicates whether the read passed Illumina's chastity filter.
	 * @return Y if failed chastity filter, N if passed
	 */
	@Override
	public char chastityCode(){return lp.parseChar(whitespaceIndex+2, 0);}

	/**
	 * Returns control bits indicating various read flags.
	 * Typically 0 for most reads, higher values indicate special conditions.
	 * @return Integer representing control bit flags
	 */
	@Override
	public int controlBits(){return lp.parseInt(whitespaceIndex+3);}

	/**
	 * Returns the barcode sequence if present in the header.
	 * May contain single index or dual indices separated by '+'.
	 * Requires a valid comment boundary; without whitespace this can read a coordinate field.
	 * @return Barcode string, or null if its term is absent after a valid comment boundary
	 */
	@Override
	public String barcode(){
		return lp.terms()<=whitespaceIndex+4 ? null : lp.parseString(whitespaceIndex+4);
	}

	/**
	 * Returns additional index information from extended header format.
	 * Available when whitespace separator is at position 7.
	 * @return Third index string or null if not available
	 */
	@Override
	public String index3(){
		return whitespaceIndex<7 ? null : lp.parseString(7);
	}

	/**
	 * Returns the position of the whitespace separator in the header.
	 * Used to determine where coordinate information ends and metadata begins.
	 * @return Term index terminated by the first space/tab, or -1 if absent
	 */
	@Override
	public int whitespaceIndex(){
		return whitespaceIndex;
	}

	/**
	 * Returns any additional information beyond standard header fields.
	 * Contains platform-specific or extended metadata.
	 * Requires a valid comment boundary; without whitespace this can read a coordinate field.
	 * @return Extra information string, or null if its term is absent after a valid comment boundary
	 */
	@Override
	public String extra(){
		return lp.terms()<=whitespaceIndex+5 ? null : lp.parseString(whitespaceIndex+5);
	}

	/** Appends the selected field without its delimiter to a caller-owned builder. */
	public ByteBuilder appendTerm(ByteBuilder bb, int term){
		return lp.appendTerm(bb, term);
	}

	/** Appends lane:tile:x:y without the instrument, run, flowcell or comment. */
	public ByteBuilder appendCoordinates(ByteBuilder bb){
		return bb.append(lane()).colon().append(tile()).colon()
		.append(xPos()).colon().append(yPos());
	}

	/** Packs lane/tile/x/y into 7/17/20/20 bits; callers must provide nonnegative values fitting those widths. */
	public long encodeCoordinates(){
		//Out-of-range fields overlap adjacent bits; this bounded-input encoding does not validate them.
		long x=lane();
		x=(x<<17)^tile();
		x=(x<<20)^xPos();
		x=(x<<20)^yPos();
		return x;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Private Fields        ----------------*/
	/*--------------------------------------------------------------*/

	private final LineParserS3 lp=new LineParserS3(':');
	/** Exposes mutable parser scratch; external mutation can invalidate this header's accessors. */
	public LineParser lp(){return lp;}
	int whitespaceIndex=-1;

}
