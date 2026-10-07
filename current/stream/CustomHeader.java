package stream;

import dna.Data;
import dna.Gene;
import parse.LineParserS2;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Stores and converts BBTools synthetic-read mapping metadata.
 * A body has nine underscore-separated fields; a complete identifier prefixes SYN_
 * and optionally joins two bodies with an ampersand and adds a mate suffix.
 * Coordinates are inclusive; start/stop are scaffold-relative or the zero-origin
 * placeholders documented by CustomHeader(Read). bbstart/bbchrom identify
 * the BBMap index, and bbstop is derived from the two coordinate systems.
 * Reference-name escaping is an ASCII format, not a general text encoding.
 * Mutable public fields and retained match arrays belong to the caller. Construction
 * from reads depends on Data's loaded genome; suffixes depend on FASTQ's global settings.
 * The older customID_old representation is separate from the SYN_ representation.
 * @author Brian Bushnell
 */
public class CustomHeader{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses a complete SYN_ identifier, selecting the mate inferred by getPairnum.
	 * @param original Nonnull identifier with a supported one-based mate suffix
	 */
	public CustomHeader(final String original){
		this(original, getPairnum(original));
	}
	
	/** Parses one body from a complete SYN_ identifier; a single body is used for either mate.
	 * The first space terminates the body. A dot match becomes null, while a dot
	 * reference name remains the literal string ".". Caught Exceptions disable
	 * FASTQ.PARSE_CUSTOM globally and may leave partial fields; this is not a validator.
	 * A local cursor checks the nine-field count before assignments; trailing empty
	 * underscore fields are omitted from the count, matching the legacy split behavior.
	 * @param original Nonnull SYN_ identifier without an initial FASTA/FASTQ marker
	 * @param rnum_ Zero-based mate, 0 or 1; explicit selection needs no suffix
	 */
	public CustomHeader(final String original, int rnum_){
		assert(rnum_==0 || rnum_==1);
		rnum=rnum_;
		
		try{
			assert(original.startsWith("SYN"));
			final int middle=original.indexOf(MIDDLE);
			final int space=original.indexOf(' ');
			//FIXED writer [stream/CustomHeader#001]: SYN_ suffixes now use modern Illumina
			//space-1:/space-2: markers (Brian, 2026-09-30), avoiding slash suffixes in rname.
			//Do not strip a legacy /1 or /2 here: it can be part of a real scaffold name.
			final String line;
			
			if(middle<0){
				line=original.substring(4, space<0 ? original.length() : space);
			}else{
				if(rnum==0){
					line=original.substring(4, middle);
				}else{
					line=original.substring(middle+1, space<0 ? original.length() : space);
				}
			}
			
			final LineParserS2 lp=new LineParserS2('_').set(line);
			int terms=(line.isEmpty() ? 1 : 0);
			for(int field=0; lp.hasMore(); field++){
				if(lp.advance()>0){terms=field+1;}
			}
			//id_start_stop_insert_strand_bbstart_bbchrom_match_encodedRname
			assert(terms==9) : "CustomHeader bodies require nine ordered metadata fields before assignment; found "+terms+": "+line;
			lp.reset();
			
			id=Long.parseLong(lp.parseString());
			start=Integer.parseInt(lp.parseString());
			stop=Integer.parseInt(lp.parseString());
			insert=Integer.parseInt(lp.parseString());
			strand=Gene.toStrand(lp.parseString());
			bbstart=Integer.parseInt(lp.parseString());
			bbchrom=Integer.parseInt(lp.parseString());
			final String matchString=lp.parseString();
			match=(matchString.equals(".") ? null : matchString.getBytes());
			rname=decodeRname(lp.parseString());
		}catch(Exception e){
			FASTQ.PARSE_CUSTOM=false;
			if(FASTQ.PARSE_CUSTOM_WARNING){
				e.printStackTrace();
				System.err.println("Turned off PARSE_CUSTOM for input "+original);
			}
		}
	}
	
	/** Stores supplied metadata, retaining match_ by reference. Three legacy arguments are ignored.
	 * @param rnum_ Zero-based mate, 0 or 1
	 * @param id_ Numeric read identity
	 * @param start_ Inclusive scaffold-relative start
	 * @param stop_ Inclusive scaffold-relative stop
	 * @param insert_ Insert-size metadata
	 * @param strand_ Shared.strandCodes2 index
	 * @param cigar_ Legacy argument, not stored or serialized
	 * @param bbstart_ Inclusive BBMap index start
	 * @param bbstop_ Legacy argument, ignored; bbstop() derives the stop from other fields
	 * @param bbchrom_ BBMap index chromosome
	 * @param bbscaffold_ Legacy argument, not stored or serialized
	 * @param match_ Match symbols without body delimiters, or null; retained without copying
	 * @param rname_ Unescaped ASCII scaffold name, or null
	 */
	public CustomHeader(int rnum_, long id_, int start_, int stop_, int insert_, int strand_, 
			String cigar_, int bbstart_, int bbstop_, int bbchrom_, int bbscaffold_, byte[] match_, String rname_){
		assert(rnum_==0 || rnum_==1);
		rnum=rnum_;
		id=id_;
		start=start_;
		stop=stop_;
		insert=insert_;
		strand=strand_;
		bbstart=bbstart_;
		bbchrom=bbchrom_;
		match=match_;
		rname=rname_;
	}
	
	/** Copies read metadata, retaining its match array and deriving scaffold-relative coordinates.
	 * With a loaded genome, the index midpoint selects a scaffold through Data.scaffoldIndex;
	 * endpoints are translated without clamping. Without a genome, start is zero and stop
	 * preserves the index span relative to that placeholder origin; rname stays null.
	 * These placeholders do not identify actual scaffold coordinates. Insert size is zero
	 * without a mate, otherwise insertSizeMapped(false).
	 * @param r Nonnull read with valid index coordinates/metadata when a genome is loaded
	 */
	public CustomHeader(Read r){
		
		rnum=r.pairnum();
		id=r.numericID;
		bbchrom=r.chrom;
		bbstart=r.start;
		int bbstop=r.stop;
		match=r.match;
		strand=r.strand();

		//[stream/CustomHeader#002 FIXED] Without a genome, retain a zero-origin span so
		//bbstop() reconstructs r.stop when FASTQ regenerates a custom name.
		start=0;
		stop=0;
		rname=null;
		int reflen=2000000000;
		
		if(Data.GENOME_BUILD>=0){
			int bbscaffold=Data.scaffoldIndex(bbchrom, (bbstart+bbstop)/2);
			byte[] name1=Data.scaffoldNames[bbchrom][bbscaffold];
			start=Data.scaffoldRelativeLoc(bbchrom, bbstart, bbscaffold);
			stop=start-bbstart+bbstop;
			rname=new String(name1);
			reflen=Data.scaffoldLengths[bbchrom][bbscaffold];
		}else{
			stop=bbstop-bbstart;
		}
		insert=(r.mate==null ? 0 : r.insertSizeMapped(false));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds this nine-field body without SYN_, a mate suffix or a line terminator.
	 * @return Newly built body; not a complete identifier for the parsing constructors */
	@Override
	public String toString(){
		ByteBuilder bb=new ByteBuilder(64);
		appendTo(bb);
		return bb.toString();
	}
	
	/** Appends id/start/stop/insert/strand/bbstart/bbchrom/match/encoded-rname in that order.
	 * Fields are underscore-separated. Null match and rname use a dot placeholder.
	 * This does not add SYN_, a mate separator/suffix or a line terminator.
	 * @param bb Nonnull caller-owned builder; existing contents are retained
	 * @return The same builder
	 */
	public ByteBuilder appendTo(ByteBuilder bb){
		bb.append(id);
		bb.append('_').append(start);
		bb.append('_').append(stop);
		bb.append('_').append(insert);
		bb.append('_').append(Shared.strandCodes2[strand]);
		bb.append('_').append(bbstart);
		bb.append('_').append(bbchrom);
		bb.append('_').append(match==null ? "." : new String(match));
		bb.append('_').append(rname==null ? "." : encodeRname(rname));
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds the older non-SYN_ custom name according to FASTQ tag/suffix globals.
	 * Returns r.id unchanged when TAG_CUSTOM is false. Otherwise starts with r.id or
	 * numericID, adds mapping fields when chrom/stop are nonnegative, and optionally
	 * translates coordinates using the loaded genome. TAG_CUSTOM_SIMPLE selects the
	 * short representation. Reference names in this legacy form are not escaped.
	 * DELETE_OLD_NAME is unsupported under assertions and otherwise clears r.id.
	 * @param r Nonnull read; valid scaffold metadata is required for genome translation
	 * @return Existing identifier or newly formatted legacy identifier
	 */
	public static String customID_old(Read r){
		if(!FASTQ.TAG_CUSTOM){return r.id;}
		
		if(FASTQ.DELETE_OLD_NAME){
			assert(false) : "Seems odd so I added this assertion.  I don't see anywhere it was being used. Use -da flag to override.";
			r.id=null;
		}

		ByteBuilder sb=new ByteBuilder();
		if(r.id==null){
			sb.append(r.numericID);
		}else{
			sb.append(r.id);
		}
		if(r.chrom>-1 && r.stop>-1){
			if(FASTQ.TAG_CUSTOM_SIMPLE){
				sb.append('_');
				sb.append(r.strand()==0 ? '+' : '-');
			}else{
				sb.append("_chr");
				sb.append(r.chrom);
				sb.append('_');
				sb.append((int)r.strand());
				sb.append('_');
				sb.append(r.start);
				sb.append('_');
				sb.append(r.stop);
			}

			if(Data.GENOME_BUILD>=0){
				final int chrom1=r.chrom;
				final int start1=r.start;
				final int stop1=r.stop;
				final int idx1=Data.scaffoldIndex(chrom1, (start1+stop1)/2);
				final byte[] name1=Data.scaffoldNames[chrom1][idx1];
				final int a1=Data.scaffoldRelativeLoc(chrom1, start1, idx1);
				final int b1=a1-start1+stop1;
				sb.append('_');
				sb.append(a1);
				if(FASTQ.TAG_CUSTOM_SIMPLE){
					sb.append('_');
					sb.append(b1);
				}
				sb.append('_');
				sb.append(new String(name1));
			}
		}
		
		if(FASTQ.ADD_PAIRNUM_TO_CUSTOM_ID){
			sb.append(' ');
			sb.append(r.pairnum()+1);
			sb.append(':');
		}else if(FASTQ.ADD_SLASH_PAIRNUM_TO_CUSTOM_ID){
			if(FASTQ.SPACE_SLASH){sb.append(' ');}
			sb.append('/');
			sb.append(r.pairnum()+1);
		}
		
		return sb.toString();
	}
	
	/** Builds a complete SYN_ name with first-mate body before second-mate body.
	 * A read with pairnum 1 must have its first mate available. A pairnum-0 read may
	 * be unpaired. Both bodies use the CustomHeader(Read) constructor's genome/array conventions.
	 * @param r Nonnull read whose pair number chooses the emitted suffix
	 * @return Complete identifier without a FASTA/FASTQ marker or line terminator
	 */
	public static String toString(Read r){
		CustomHeader h1=new CustomHeader(r);
		CustomHeader h2=(r.mate==null ? null : new CustomHeader(r.mate));
		final String s;
		if(r.pairnum()==0){
			s=toString(h1, h2, r.pairnum());
		}else{
			s=toString(h2, h1, r.pairnum());
		}
		return s;
	}
	
	/** Builds SYN_ plus one or two bodies in the supplied order, separated by an ampersand.
	 * Either ADD_PAIRNUM_TO_CUSTOM_ID or ADD_SLASH_PAIRNUM_TO_CUSTOM_ID enables the
	 * modern Illumina suffix: space, one-based mate, colon. SPACE_SLASH does not alter
	 * SYN_ output. With both suffix flags false no mate suffix is appended. Body rnum
	 * fields do not determine ordering or the suffix.
	 * @param h1 Nonnull first body
	 * @param h2 Second body, or null for a single-body identifier
	 * @param rnum Zero-based mate, 0 or 1, used only for the suffix
	 * @return New complete identifier without a line terminator
	 */
	public static String toString(CustomHeader h1, CustomHeader h2, int rnum){
		ByteBuilder bb=new ByteBuilder();
		bb.append("SYN_");
		h1.appendTo(bb);
		if(h2!=null){
			bb.append('&');
			h2.appendTo(bb);
		}
		
		if(FASTQ.ADD_PAIRNUM_TO_CUSTOM_ID || FASTQ.ADD_SLASH_PAIRNUM_TO_CUSTOM_ID){
			bb.append(' ');
			bb.append(rnum+1);
			bb.append(':');
		}
		return bb.toString();
	}
	
	/** Decodes an ASCII reference-name token: dollar to space, left brace to underscore,
	 * and !a/!b/!d/!! to ampersand/left brace/dollar/exclamation mark.
	 * A literal dot stays a dot. The caller supplies a complete token produced by encodeRname.
	 * @param rname Nonnull encoded ASCII name
	 * @return Decoded name in a new string
	 */
	public static String decodeRname(String rname){
		ByteBuilder bb=new ByteBuilder(rname.length());
		char prev='.';
		for(int i=0; i<rname.length(); i++){
			char c=rname.charAt(i);
			if(prev=='!'){
				bb.append(escapeArray[c]);
				prev='.';
			}else if(c!='!'){
				bb.append(decodeArray[c]);
				prev=c;
			}else{
				prev=c;
			}
		}
		return bb.toString();
	}
	
	/** Escapes the six reserved ASCII characters used by the SYN_ name format.
	 * Spaces become dollars, underscores become left braces, and !/ampersand/dollar/
	 * left brace become !!/!a/!d/!b respectively. Other characters are passed through;
	 * this does not sanitize tabs/newlines or provide a Unicode encoding.
	 * @param rname Nonnull unescaped ASCII name without line/tab delimiters
	 * @return Original string if no escaping is needed, otherwise a new encoded string
	 */
	public static String encodeRname(String rname){
		
		int found=0;
		for(int i=0; i<rname.length(); i++){
			char c=rname.charAt(i);
			if(c=='!' || c=='&' || c==' ' || c=='$' || c=='{' || c=='_'){found++;}
		}
		if(found<1){return rname;}
		
		ByteBuilder bb=new ByteBuilder(rname.length()+found);
		for(int i=0; i<rname.length(); i++){
			char c=rname.charAt(i);
			if(c=='!'){bb.append("!!");}
			else if(c==' '){bb.append("$");}
			else if(c=='_'){bb.append("{");}
			else if(c=='&'){bb.append("!a");}
			else if(c=='$'){bb.append("!d");}
			else if(c=='{'){bb.append("!b");}
			else{bb.append(c);}
		}
		return bb.toString();
	}
	
	/** Infers a zero-based mate from the character following the last space or slash.
	 * A space takes precedence and may be followed by a slash. Supported suffixes
	 * include space-1:, space-/1 and /1, with 2 for the second mate. This helper does
	 * not validate the remaining suffix and is unsuitable for suffix-free names.
	 * @param s Nonnull identifier with a supported one-based mate suffix
	 * @return First suffix digit minus '1', normally 0 or 1
	 */
	public static int getPairnum(String s){
		int idx=s.lastIndexOf(' ');
		if(idx<0){idx=s.lastIndexOf('/');}
		assert(idx>0) : s;
		idx++;
		if(s.charAt(idx)=='/'){idx++;}
		return s.charAt(idx)-'1';
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Zero-based selected mate; not included in the body. */
	public int rnum;
	/** Numeric read identity. */
	public long id;
	/** Inclusive scaffold-relative start, or zero when no genome was available. */
	public int start;
	/** Inclusive scaffold-relative stop, or the end of the zero-origin placeholder span. */
	public int stop;
	/** Insert-size metadata; Read construction delegates to its mapped-insert helper. */
	public int insert;
	/** Strand code indexing Shared.strandCodes2. */
	public int strand;
	/** Inclusive BBMap index start. */
	public int bbstart;
	/** Derives inclusive BBMap index stop by preserving the stored scaffold-relative span. */
	public int bbstop(){return bbstart-start+stop;}
	/** BBMap index chromosome. */
	public int bbchrom;
	/** Match symbols, or null; constructors from a Read or supplied array retain its reference. */
	public byte[] match;
	/** Unescaped scaffold name; null writes a dot, which parses back as a literal dot. */
	public String rname;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Historical escape marker; codec branches currently use literal characters. */
	private static final char ESCAPE='!';
	/** Separator between the two mate bodies. */
	private static final char MIDDLE='&';

	/** ASCII byte mapping outside an escape sequence. */
	private static final byte[] decodeArray;
	/** ASCII byte mapping immediately after an escape marker. */
	private static final byte[] escapeArray;
	
	static{
		decodeArray=new byte[128];
		escapeArray=new byte[128];
		for(int i=0; i<128; i++){
			decodeArray[i]=escapeArray[i]=(byte)i;
			if(i=='$'){decodeArray[i]=' ';}
			else if(i=='{'){decodeArray[i]='_';}
			else if(i=='a'){escapeArray[i]='&';}
			else if(i=='b'){escapeArray[i]='{';}
			else if(i=='d'){escapeArray[i]='$';}
		}
	}
	
}
