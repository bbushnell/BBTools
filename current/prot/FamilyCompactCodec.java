package prot;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;

import parse.LineParser1;
import structures.ByteBuilder;

/** Encodes and decodes ten-line compact family blocks without file I/O.
 * The header has one magic/version line followed by metadata and blank padding.
 * Family blocks use fam, gate, rank, gc, seq, ddl, mono, pos, res1 and res2 order.
 * The last three tags have no payload in version 1; empty families have only
 * their fam payload and nine bare tags. DDL registers use padded three-byte
 * uint16 A48; floating-point fields use variable-length raw-bit A48 instead.
 * Header parsing validates syntax, not required package metadata or resource
 * bindings. Whole-file identity/census checks belong to FamilyCompactLoader;
 * the callable package validates its geometry, matrix and resource policy.
 * Input lines are borrowed for the call. Parsers and output builders are
 * caller-owned mutable scratch and must not be shared by simultaneous calls.
 * Append methods retain any existing output prefix and do not roll it back on
 * failure. Parsed payload arrays are newly allocated and handed to the sink.
 * @author Yelan
 */
public final class FamilyCompactCodec{

	/** Static wire-format operations only. */
	private FamilyCompactCodec(){}

	/** Physical block and ordinary input-batch sizes; the header is itself one block. */
	public static final int LINES_PER_RECORD=10,LINES_PER_BATCH=200,FAMILIES_PER_BATCH=20;
	/** Fixed uint16 register payload and stored monomer-vector dimensions. */
	public static final int DDL_WIDTH=512,DDL_BYTES=3*DDL_WIDTH,MONOMER_WIDTH=20;
	/** Exact first-line magic and supported format version. */
	public static final String MAGIC="#bbgap_family_compact",VERSION="1";

	/** Explicit alignment-calibration availability, independent of scalar-gate bits. */
	public enum AlignmentCalibration{
		ENABLED("alignment_enabled"),
		UNAVAILABLE("alignment_unavailable");
		/** Stores the exact wire token. */
		AlignmentCalibration(String label_){label=label_;}
		public final String label;
	}

	/** Detached, unmodifiable metadata map preserving the input map's iteration order.
	 * Values may contain '='; keys may not. Terms must be nonempty and contain no
	 * tab, carriage return or newline. Required keys and semantic values are checked
	 * by the consuming loader/package, not this container.
	 */
	public static final class Header{
		/** Copies and checks lexical field constraints without retaining the supplied map. */
		public Header(Map<String,String> fields_){
			require(fields_!=null && !fields_.isEmpty(),"Compact header requires metadata fields");
			final LinkedHashMap<String,String> copy=new LinkedHashMap<String,String>();
			for(Map.Entry<String,String> entry : fields_.entrySet()){
				final String key=entry.getKey(),value=entry.getValue();
				require(key!=null && value!=null && key.indexOf('=')<0 && cleanTerm(key) && cleanTerm(value),
					"Compact header fields must be nonempty key=value tokens");
				require(copy.put(key,value)==null,"Duplicate compact header key");
			}
			fields=Collections.unmodifiableMap(copy);
		}
		/** Returns the stored value, or null when the key is absent. */
		public String get(String key){return fields.get(key);}
		public final Map<String,String> fields;
	}

	/** Shallow value carrier for one physical family, including zero-member families.
	 * Array fields are borrowed directly and remain publicly mutable. Construction
	 * does not validate the values; appendRecord checks them before or during output.
	 * Keep arrays stable while writing. Rank, stable ID and physical file position
	 * are distinct identities; this class does not enforce whole-file uniqueness.
	 */
	public static final class Record{
		/** Stores metadata and array references without copying, normalization or validation. */
		public Record(String familyId_,int rank_,int stableId_,long members_,long residues_,long validWords_,
				float cosineLower_,float identityLower_,float blosumLower_,float cosineAccept_,
				float identityUpper_,float blosumFinal_,int scalarMask_,
				AlignmentCalibration alignmentCalibration_,double[] scalarBounds_,
				double lengthMedian_,double kdMean_,long measured_,double meanGc_,double sampleSd_,
				byte[] consensusRaw_,char[] ddl_,float[] monomers_){
			familyId=familyId_; rank=rank_; stableId=stableId_; members=members_; residues=residues_;
			validWords=validWords_; cosineLower=cosineLower_; identityLower=identityLower_;
			blosumLower=blosumLower_; cosineAccept=cosineAccept_; identityUpper=identityUpper_;
			blosumFinal=blosumFinal_; scalarMask=scalarMask_;
			alignmentCalibration=alignmentCalibration_; scalarBounds=scalarBounds_;
			lengthMedian=lengthMedian_; kdMean=kdMean_; measured=measured_; meanGc=meanGc_;
			sampleSd=sampleSd_; consensusRaw=consensusRaw_; ddl=ddl_; monomers=monomers_;
		}
		public final String familyId;
		public final int rank,stableId,scalarMask;
		public final long members,residues,validWords,measured;
		public final float cosineLower,identityLower,blosumLower,cosineAccept,identityUpper,blosumFinal;
		public final AlignmentCalibration alignmentCalibration;
		public final double[] scalarBounds;
		public final double lengthMedian,kdMean,meanGc,sampleSd;
		public final byte[] consensusRaw;
		public final char[] ddl;
		public final float[] monomers;
	}

	/** Receives exactly one callback after a complete record passes codec validation.
	 * Live arrays are newly decoded, independent of the input lines, and may be
	 * retained by the sink. No intermediate Record wrapper is needed by custom sinks.
	 * Sink exceptions propagate; the codec does not undo sink-side effects.
	 * @param <T> Caller-selected record representation
	 */
	public interface RecordSink<T>{
		/** Builds a zero-member record from its retained identities. */
		T empty(String familyId,int rank,int stableId);
		/** Takes ownership of decoded live payload arrays and scalar metadata. */
		T live(String familyId,int rank,int stableId,long members,long residues,long validWords,
				float cosineLower,float identityLower,float blosumLower,float cosineAccept,
				float identityUpper,float blosumFinal,int scalarMask,
				AlignmentCalibration alignmentCalibration,double[] scalarBounds,
				double lengthMedian,double kdMean,long measured,double meanGc,double sampleSd,
				byte[] consensusRaw,char[] ddl,float[] monomers);
	}

	/** Header line zero is fixed magic/version; lines one through nine are tab-packed key=value terms.
	 * Blank physical lines are right padding only, preserving the ten-line boundary.
	 * @param lines Nonnull physical lines without line terminators; no rows are retained.
	 * @param offset First header line within lines.
	 * @param lp Caller-owned tab-delimited parser, reused and left at its last parsed row.
	 * @return Detached syntactically validated metadata; required keys are not inferred.
	 */
	public static Header parseHeader(byte[][] lines,int offset,LineParser1 lp){
		require(lines!=null && lp!=null && offset>=0 && offset<=lines.length-LINES_PER_RECORD,
			"Incomplete ten-line compact header");
		for(int i=0; i<LINES_PER_RECORD; i++){require(lines[offset+i]!=null,"Null compact header line");}
		lp.set(lines[offset]); require(lp.terms()==2 && lp.termEquals(MAGIC,0) && lp.termEquals(VERSION,1),
			"Expected compact family magic/version");
		final LinkedHashMap<String,String> fields=new LinkedHashMap<String,String>();
		boolean padded=false;
		for(int line=1; line<LINES_PER_RECORD; line++){
			final byte[] row=lines[offset+line];
			if(row.length==0){padded=true; continue;}
			require(!padded,"Compact header metadata cannot follow a padded blank line");
			lp.set(row);
			for(int term=0; term<lp.terms(); term++){
				final String field=lp.parseString(term);
				final int equals=field.indexOf('=');
				require(equals>0 && equals<field.length()-1,"Compact header fields must be key=value");
				final String key=field.substring(0,equals);
				require(fields.put(key,field.substring(equals+1))==null,"Duplicate compact header key");
			}
		}
		return new Header(fields);
	}

	/** Appends magic/version, one packed metadata line and eight blank padding lines.
	 * Preserves field iteration order and the existing output prefix; returns out.
	 */
	public static ByteBuilder appendHeader(Header h,ByteBuilder out){
		require(h!=null && out!=null,"Missing compact header/output");
		out.append(MAGIC).tab().append(VERSION).nl();
		boolean first=true;
		for(Map.Entry<String,String> entry : h.fields.entrySet()){
			if(!first){out.tab();}
			out.append(entry.getKey()).append('=').append(entry.getValue());
			first=false;
		}
		for(int i=1; i<LINES_PER_RECORD; i++){out.nl();}
		return out;
	}

	/** Decodes into a new Record with owned payload arrays, including inert arrays for empties. */
	public static Record parseRecord(byte[][] lines,int offset,LineParser1 lp){
		return parseRecord(lines,offset,lp,RECORD_SINK);
	}

	/** Validates one complete ten-line family block and returns the sink's result.
	 * Checks local tags, dimensions, encodings and numeric constraints, not header
	 * policy or uniqueness across records. Live payload arrays are allocated once
	 * for the sink; the codec retains neither the sink nor the input lines.
	 * @param <T> Sink result type
	 * @param lines Physical lines without terminators.
	 * @param offset First fam line within lines.
	 * @param lp Exclusively owned tab-delimited parser; its state is overwritten.
	 * @param sink Consumer of validated empty/live records.
	 * @return Exactly the object returned by the selected sink callback.
	 */
	public static <T> T parseRecord(byte[][] lines,int offset,LineParser1 lp,RecordSink<T> sink){
		require(lines!=null && lp!=null && offset>=0 && offset<=lines.length-LINES_PER_RECORD,
			"Incomplete ten-line compact family record");
		require(sink!=null,"Missing compact record sink");
		for(int i=0; i<LINES_PER_RECORD; i++){require(lines[offset+i]!=null,"Null compact family line");}
		lp.set(lines[offset]); require(lp.terms()==7 && lp.termEquals("fam",0),"Expected compact fam row");
		final int rank=integer(lp,1);
		final String familyId=lp.parseString(2);
		final int stableId=integer(lp,3);
		final long members=nonnegativeLong(lp,4),residues=nonnegativeLong(lp,5),validWords=nonnegativeLong(lp,6);
		require(!familyId.isEmpty() && familyId.indexOf('\n')<0 && familyId.indexOf('\r')<0,
			"Compact family ID must be nonempty");
		require(validWords<=residues,"Compact family valid words exceed residues");

		if(members==0){
			require(residues==0 && validWords==0,"Empty compact families require zero counts");
			reserved(lines[offset+1],"gate");
			reserved(lines[offset+2],"rank");
			reserved(lines[offset+3],"gc");
			reserved(lines[offset+4],"seq");
			reserved(lines[offset+5],"ddl");
			reserved(lines[offset+6],"mono");
			reserved(lines[offset+7],"pos");
			reserved(lines[offset+8],"res1");
			reserved(lines[offset+9],"res2");
			return sink.empty(familyId,rank,stableId);
		}

		lp.set(lines[offset+1]); require(lp.terms()==19 && lp.termEquals("gate",0),"Expected compact gate row");
		final float cosineLower=FamilyGateCodec.readLower(lp,1),identityLower=FamilyGateCodec.readLower(lp,2),
			blosumLower=FamilyGateCodec.readLower(lp,3),cosineAccept=FamilyGateCodec.readCosineAccept(lp,4),
			identityUpper=FamilyGateCodec.readUpper(lp,5),blosumFinal=FamilyGateCodec.readFinalBlosum(lp,6);
		final int scalarMask=integer(lp,7);
		final AlignmentCalibration alignment=alignmentCalibration(lp,8);
		require(scalarMask>=0 && scalarMask<32,"Compact scalar mask must be five bits");
		final double[] scalarBounds=new double[10];
		for(int i=0; i<scalarBounds.length; i++){scalarBounds[i]=readDoubleA48(lp,9+i);}
		for(int i=0; i<5; i++){
			final double lo=scalarBounds[2*i],hi=scalarBounds[2*i+1];
			require((scalarMask&(1<<i))!=0 ? members>=30 && lo<hi : lo==0 && hi==0,
				"Compact scalar gates need N>=30 and increasing bounds or exact inert zeros");
		}

		lp.set(lines[offset+2]); require(lp.terms()==3 && lp.termEquals("rank",0),
			"Expected compact rank row");
		final double lengthMedian=readDoubleA48(lp,1),kdMean=readDoubleA48(lp,2);
		require(lengthMedian>=40 && lengthMedian<=200 && lengthMedian*2==Math.rint(lengthMedian*2),
			"Compact length median must be a finite integer or half-integer in 40..200");

		lp.set(lines[offset+3]); require(lp.terms()==4 && lp.termEquals("gc",0),"Expected compact gc row");
		final long measured=nonnegativeLong(lp,1);
		require(measured==members,"Compact GC measurement count must match fam members");
		final double meanGc=fraction(lp,2),sampleSd=fraction(lp,3);
		require(members>1 || sampleSd==0,"Singleton and empty compact GC sample SD must be zero");

		lp.set(lines[offset+4]); require(lp.terms()==2 && lp.termEquals("seq",0),"Expected compact seq row");
		final byte[] consensus=termBytes(lp,1);
		require(consensus.length>0,"Compact consensus must be nonempty");
		require(supportedConsensus(consensus),"Compact consensus contains unsupported amino acid bytes");

		final char[] ddl=new char[DDL_WIDTH];
		parseFixedDdl(lines[offset+5],ddl);

		lp.set(lines[offset+6]); require(lp.terms()==MONOMER_WIDTH+1 && lp.termEquals("mono",0),
			"Expected compact mono row");
		final float[] monomers=new float[MONOMER_WIDTH];
		for(int i=0; i<monomers.length; i++){
			monomers[i]=readFloatA48(lp,i+1);
			require(Float.isFinite(monomers[i]) && monomers[i]>=0 && monomers[i]<=1,
				"Compact monomer frequencies must be finite fractions");
		}
		reserved(lines[offset+7],"pos");
		reserved(lines[offset+8],"res1");
		reserved(lines[offset+9],"res2");
		return sink.live(familyId,rank,stableId,members,residues,validWords,cosineLower,
			identityLower,blosumLower,cosineAccept,identityUpper,blosumFinal,scalarMask,
			alignment,scalarBounds,lengthMedian,kdMean,measured,meanGc,sampleSd,consensus,ddl,monomers);
	}

	/** Appends one validated family block and its final newline, returning out.
	 * Empty records serialize only identity/zero counts and bare payload tags.
	 * Gate checks occur during writing, so a rejected record can leave a partial
	 * prefix in out. Callers discard or clear failed output before reuse.
	 */
	public static ByteBuilder appendRecord(Record r,ByteBuilder out){
		require(r!=null && out!=null,"Missing compact record/output");
		validateRecord(r);
		out.append("fam").tab().append(r.rank).tab().append(r.familyId).tab().append(r.stableId)
			.tab().append(r.members).tab().append(r.residues).tab().append(r.validWords).nl();
		if(r.members==0){return out.append("gate\nrank\ngc\nseq\nddl\nmono\npos\nres1\nres2\n");}
		out.append("gate").tab();
		FamilyGateCodec.appendLower(out,r.cosineLower).tab();
		FamilyGateCodec.appendLower(out,r.identityLower).tab();
		FamilyGateCodec.appendLower(out,r.blosumLower).tab();
		FamilyGateCodec.appendCosineAccept(out,r.cosineAccept).tab();
		FamilyGateCodec.appendUpper(out,r.identityUpper).tab();
		FamilyGateCodec.appendFinalBlosum(out,r.blosumFinal).tab().append(r.scalarMask).tab()
			.append(r.alignmentCalibration.label);
		for(double bound : r.scalarBounds){out.tab().appendDoubleA48(bound);}
		out.nl().append("rank").tab();
		out.appendDoubleA48(r.lengthMedian).tab().appendDoubleA48(r.kdMean);
		out.nl()
			.append("gc").tab().append(r.measured).tab().appendDoubleA48(r.meanGc).tab()
			.appendDoubleA48(r.sampleSd).nl()
			.append("seq").tab().append(r.consensusRaw).nl();
		appendFixedDdl(r.ddl,out).nl().append("mono");
		for(float value : r.monomers){out.tab().appendFloatA48(value);}
		return out.append("\npos\nres1\nres2\n");
	}

	/** Checks record metadata/payloads; finite gate checks are delegated to gate writers.
	 * Zero-member records need only consistent fam fields; other payloads are not emitted.
	 */
	private static void validateRecord(Record r){
		require(r.familyId!=null && cleanTerm(r.familyId),"Compact family ID must be nonempty");
		require(r.rank>=0 && r.stableId>=0 && r.members>=0 && r.residues>=0 &&
			r.validWords>=0 && r.validWords<=r.residues,"Compact fam counts must be consistent");
		if(r.members==0){
			require(r.residues==0 && r.validWords==0,"Empty compact families require zero counts");
			return;
		}
		require(r.alignmentCalibration!=null,"Compact alignment calibration is missing");
		require(r.scalarBounds!=null && r.scalarBounds.length==10,"Compact scalar bounds require ten values");
		require(r.scalarMask>=0 && r.scalarMask<32,"Compact scalar mask must be five bits");
		for(int i=0; i<5; i++){
			final double lo=r.scalarBounds[2*i],hi=r.scalarBounds[2*i+1];
			require(Double.isFinite(lo) && Double.isFinite(hi),
				"Compact scalar bounds must be finite");
			require((r.scalarMask&(1<<i))!=0 ? r.members>=30 && lo<hi : lo==0 && hi==0,
				"Compact scalar gates need N>=30 and increasing bounds or exact inert zeros");
		}
		require(Double.isFinite(r.lengthMedian) && r.lengthMedian>=40 && r.lengthMedian<=200 &&
			r.lengthMedian*2==Math.rint(r.lengthMedian*2),
			"Compact length median must be a finite integer or half-integer in 40..200");
		require(Double.isFinite(r.kdMean),"Compact KD mean must be finite");
		require(r.measured==r.members,"Compact GC measurement count must match fam members");
		require(Double.isFinite(r.meanGc) && r.meanGc>=0 && r.meanGc<=1 &&
			Double.isFinite(r.sampleSd) && r.sampleSd>=0 && r.sampleSd<=1,
			"Compact GC values must be finite fractions");
		require(r.members>1 || r.sampleSd==0,"Singleton compact GC sample SD must be zero");
		require(r.consensusRaw!=null && r.consensusRaw.length>0 && supportedConsensus(r.consensusRaw),
			"Compact consensus must be nonempty and supported");
		require(r.ddl!=null && r.ddl.length==DDL_WIDTH,"Compact DDL requires 512 registers");
		require(r.monomers!=null && r.monomers.length==MONOMER_WIDTH,
			"Compact monomers require 20 frequencies");
		for(float value : r.monomers){
			require(Float.isFinite(value) && value>=0 && value<=1,
				"Compact monomer frequencies must be finite fractions");
		}
	}

	/** Appends the ddl tag, a tab and exactly 1536 payload bytes, without a newline.
	 * Each uint16 uses three high-to-low radix-64 digits offset by ASCII 48;
	 * zero is 000 and 65535 is ?oo. Returns the caller's existing builder.
	 */
	public static ByteBuilder appendFixedDdl(char[] ddl,ByteBuilder out){
		require(ddl!=null && ddl.length==DDL_WIDTH && out!=null,"Compact DDL requires 512 registers");
		out.append("ddl").tab();
		for(char c : ddl){appendFixedDdlValue(c,out);}
		return out;
	}

	/** Emits one uint16 as three digits, retaining leading zero digits. */
	private static void appendFixedDdlValue(char c,ByteBuilder out){
		out.append((byte)(48+(c>>12)));
		out.append((byte)(48+((c>>6)&63)));
		out.append((byte)(48+(c&63)));
	}

	/** Validates a fixed-size ddl row and fills the caller's 512-register array. */
	private static void parseFixedDdl(byte[] line,char[] ddl){
		require(line.length==4+DDL_BYTES && line[0]=='d' && line[1]=='d' && line[2]=='l' &&
			line[3]=='\t',"Compact DDL row must contain exactly 1536 payload bytes");
		int pos=4;
		for(int i=0; i<ddl.length; i++,pos+=3){
			final int a=line[pos]-48,b=line[pos+1]-48,c=line[pos+2]-48;
			require(a>=0 && a<16 && b>=0 && b<64 && c>=0 && c<64,
				"Invalid fixed-width compact DDL alphabet");
			ddl[i]=(char)((a<<12)|(b<<6)|c);
		}
	}

	/** Decodes validated unsigned32 raw bits; the caller enforces the float's numeric domain. */
	private static float readFloatA48(LineParser1 lp,int term){
		final long bits=a48(lp,term,0xffffffffL);
		return Float.intBitsToFloat((int)bits);
	}

	/** Decodes a raw-bit A48 double and rejects nonfinite values. */
	private static double readDoubleA48(LineParser1 lp,int term){
		require(validA48_64(lp,term),"Invalid compact A48 f64 field");
		final double value=lp.parseDoubleA48(term);
		require(Double.isFinite(value),"Compact A48 f64 field must be finite");
		return value;
	}

	/** Validates a minimal unsigned token of at most six digits before narrowing. */
	private static long a48(LineParser1 lp,int term,long max){
		final int length=lp.setBounds(term),start=lp.a();
		final byte[] line=lp.line();
		require(length>0 && length<=6 && (length==1 || line[start]!='0'),"Noncanonical A48 token");
		for(int i=start; i<lp.b(); i++){require(line[i]>=48 && line[i]<=111,"Invalid A48 alphabet");}
		final long value=lp.parseLongA48(term);
		require(value>=0 && value<=max,"A48 value exceeds field range");
		return value;
	}

	/** Checks the variable-width long encoding, including sign-extended leading digits.
	 * This lexical check does not establish that the decoded double is finite.
	 */
	private static boolean validA48_64(LineParser1 lp,int term){
		final int length=lp.setBounds(term);
		if(length<1 || length>11){return false;}
		final byte[] line=lp.line();
		final int first=line[lp.a()]-48;
		if(length>1 && line[lp.a()]=='0'){return false;}
		if(length==11 && first>=8 && first<56){return false;}
		for(int i=lp.a(); i<lp.b(); i++){
			if(line[i]<48 || line[i]>111){return false;}
		}
		return true;
	}

	/** Accepts only the two explicit wire labels, with no inferred default. */
	private static AlignmentCalibration alignmentCalibration(LineParser1 lp,int term){
		if(lp.termEquals(AlignmentCalibration.ENABLED.label,term)){
			return AlignmentCalibration.ENABLED;
		}
		if(lp.termEquals(AlignmentCalibration.UNAVAILABLE.label,term)){
			return AlignmentCalibration.UNAVAILABLE;
		}
		throw new IllegalArgumentException("Unknown compact alignment calibration");
	}

	/** Requires a bare tag with no tab or payload, including for empty-family lines. */
	private static void reserved(byte[] line,String tag){
		require(tagEquals(line,tag),"Compact reserved row must be the bare tag");
	}

	/** Compares the complete line with the expected ASCII tag. */
	private static boolean tagEquals(byte[] line,String tag){
		if(line.length!=tag.length()){return false;}
		for(int i=0; i<line.length; i++){
			if(line[i]!=tag.charAt(i)){return false;}
		}
		return true;
	}

	/** Copies a field so retained sequence bytes do not borrow the input line. */
	private static byte[] termBytes(LineParser1 lp,int term){
		lp.setBounds(term);
		final byte[] out=new byte[lp.b()-lp.a()];
		System.arraycopy(lp.line(),lp.a(),out,0,out.length);
		return out;
	}

	/** Accepts standard uppercase amino acids plus X, B, Z and J; no stops or lowercase. */
	private static boolean supportedConsensus(byte[] consensus){
		for(byte b : consensus){
			switch(b){
				case 'A': case 'R': case 'N': case 'D': case 'C':
				case 'E': case 'Q': case 'G': case 'H': case 'I':
				case 'L': case 'K': case 'M': case 'F': case 'P':
				case 'S': case 'T': case 'W': case 'Y': case 'V':
				case 'X': case 'B': case 'Z': case 'J':
					break;
				default:
					return false;
			}
		}
		return true;
	}

	/** Reads a finite raw-bit double in the inclusive range zero through one. */
	private static double fraction(LineParser1 lp,int term){
		final double value=readDoubleA48(lp,term);
		require(value>=0 && value<=1,"Compact fraction outside 0..1");
		return value;
	}

	/** Reads a nonnegative decimal count bounded by Integer.MAX_VALUE. */
	private static int integer(LineParser1 lp,int term){
		final long x=nonnegativeLong(lp,term);
		require(x<=Integer.MAX_VALUE,"Compact integer overflow");
		return (int)x;
	}

	/** Reads a decimal long and rejects negative values. */
	private static long nonnegativeLong(LineParser1 lp,int term){
		final long value=signedLong(lp,term);
		require(value>=0,"Compact count/rank must be nonnegative");
		return value;
	}

	/** Parses optional minus plus decimal digits, checking overflow before every operation. */
	private static long signedLong(LineParser1 lp,int term){
		final int length=lp.setBounds(term);
		require(length>0,"Empty compact decimal integer");
		final byte[] line=lp.line();
		int pos=lp.a();
		final boolean negative=line[pos]=='-';
		if(negative){pos++;}
		require(pos<lp.b(),"Missing compact decimal digits");
		long value=0;
		final long limit=negative ? Long.MIN_VALUE : -Long.MAX_VALUE;
		for(; pos<lp.b(); pos++){
			final int digit=line[pos]-'0';
			require(digit>=0 && digit<=9 && value>=limit/10,"Invalid or overflowing compact decimal integer");
			value*=10;
			require(value>=limit+digit,"Compact decimal integer overflow");
			value-=digit;
		}
		return negative ? value : -value;
	}

	/** Unconditional input validation, independent of JVM assertions. */
	private static void require(boolean condition,String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Tests nonempty delimiter-free text; does not restrict the remaining character set. */
	private static boolean cleanTerm(String s){
		return !s.isEmpty() && s.indexOf('\t')<0 && s.indexOf('\n')<0 && s.indexOf('\r')<0;
	}

	/** Default sink retains decoded arrays directly; empty records receive inert arrays. */
	private static final RecordSink<Record> RECORD_SINK=new RecordSink<Record>(){
		@Override
		public Record empty(String familyId,int rank,int stableId){
			return new Record(familyId,rank,stableId,0,0,0,0,0,0,0,0,0,0,
				AlignmentCalibration.UNAVAILABLE,new double[10],0,0,0,0,0,new byte[0],
				new char[DDL_WIDTH],new float[MONOMER_WIDTH]);
		}
		@Override
		public Record live(String familyId,int rank,int stableId,long members,long residues,
				long validWords,float cosineLower,float identityLower,float blosumLower,
				float cosineAccept,float identityUpper,float blosumFinal,int scalarMask,
				AlignmentCalibration alignmentCalibration,double[] scalarBounds,
				double lengthMedian,double kdMean,long measured,double meanGc,double sampleSd,
				byte[] consensusRaw,char[] ddl,float[] monomers){
			return new Record(familyId,rank,stableId,members,residues,validWords,cosineLower,
				identityLower,blosumLower,cosineAccept,identityUpper,blosumFinal,scalarMask,
				alignmentCalibration,scalarBounds,lengthMedian,kdMean,measured,meanGc,sampleSd,
				consensusRaw,ddl,monomers);
		}
	};
}
