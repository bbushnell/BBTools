package prot;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import parse.LineParser1;
import structures.ByteBuilder;

/** Five-line family codec with explicit dimer or projected-monomer storage. Parsers are caller-owned and thread-local.
 * @author Sayu
 */
public final class FamilyReferenceCodec{

	private FamilyReferenceCodec(){}
	public static final int LINES_PER_FAMILY=5,DDL_WIDTH=1024,SPECTRUM_WIDTH=400,MONOMER_WIDTH=20;
	public static final String SCHEMA="family_reference_v7",
		FORMULA_VERSION="k5_identity30_blosum_uncapped_v4",SPECTRUM_ENCODING="a48_f32";
	public static final String MONOMER_SCHEMA="family_reference_v8";
	/** New geometry is explicit; unchanged v7/v8 headers retain their exact historical bytes. */
	public static final String METRIC_SCHEMA="family_reference_v9",METRIC_MONOMER_SCHEMA="family_reference_v10";
	public static final String ALIGNMENT_SCHEMA="family_reference_v11",ALIGNMENT_MONOMER_SCHEMA="family_reference_v12";
	public static final String MATRIX_SCHEMA="family_reference_v13",MATRIX_MONOMER_SCHEMA="family_reference_v14";
	public static final String PAIRED_SCHEMA="family_reference_v15",PAIRED_MONOMER_SCHEMA="family_reference_v16";
	private static final byte[] MONOMER_LABEL="#spectrum\t20\ta48_f32".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] COMPACT_DDL_LABEL="#ddl\t512\tA48".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] ALIGNMENT_GATES_LABEL=("#gates\tcosine_lower\tidentity_lower\tblosum_lower\tcosine_accept"+
		"\tidentity_upper\tblosum_final\tscalar_mask\talignment_state\thydro_min\thydro_max\tcharge_min\tcharge_max"+
		"\taromatic_min\taromatic_max\tsmall_min\tsmall_max\tpi_min\tpi_max").getBytes(StandardCharsets.US_ASCII);
	private static final byte[][] HEADER_LABELS={
		("#gates\tcosine_lower\tidentity_lower\tblosum_lower\tcosine_accept\tidentity_upper"+
		"\tblosum_final\tscalar_mask\thydro_min\thydro_max\tcharge_min\tcharge_max\taromatic_min"+
		"\taromatic_max\tsmall_min\tsmall_max\tpi_min\tpi_max").getBytes(StandardCharsets.US_ASCII),
		"#consensus\tamino_acids".getBytes(StandardCharsets.US_ASCII),
		"#ddl\t1024\tA48".getBytes(StandardCharsets.US_ASCII),
		"#spectrum\t400\ta48_f32".getBytes(StandardCharsets.US_ASCII)
	};

	public static final class Header{
		public Header(String schema,int k,int width,long seed,int totalFamilies,int liveFamilies,int emptyFamilies,
				String formulaVersion,String spectrumEncoding){
			this(schema,k,width,seed,totalFamilies,liveFamilies,emptyFamilies,formulaVersion,spectrumEncoding,0);
		}
		public Header(String schema,int k,int width,long seed,int totalFamilies,int liveFamilies,int emptyFamilies,
				String formulaVersion,String spectrumEncoding,double scalarZ){
			this(schema,k,width,seed,totalFamilies,liveFamilies,emptyFamilies,formulaVersion,spectrumEncoding,
				scalarZ,FamilyAlignmentMetric.LEGACY);
		}
		/** A new schema requires an explicit pinned metric and positive span. */
		public Header(String schema,int k,int width,long seed,int totalFamilies,int liveFamilies,int emptyFamilies,
				String formulaVersion,String spectrumEncoding,double scalarZ,FamilyAlignmentMetric alignmentMetric){
			require(Double.isFinite(scalarZ) && scalarZ>=0,
				"Scalar z must be finite and nonnegative; zero disables package gates");
			this.scalarZ=scalarZ;
			final boolean bound=metricSchema(schema);
			require((SCHEMA.equals(schema) || MONOMER_SCHEMA.equals(schema) || bound) &&
				FORMULA_VERSION.equals(formulaVersion) && SPECTRUM_ENCODING.equals(spectrumEncoding),
				"Unsupported family reference schema/formula/encoding");
			require(alignmentMetric!=null && (bound ? pinnedMetric(alignmentMetric.mode) :
				alignmentMetric.same(FamilyAlignmentMetric.LEGACY)),
				"Reference schema and alignment metric disagree; new geometry needs an explicit new package");
			require(pairedSchema(schema) ? alignmentMetric.minimumPaired>0 : alignmentMetric.minimumPaired==0,
				"Paired minima require an explicit paired-bound reference schema");
			this.alignmentMetric=alignmentMetric;
			require(matrixSchema(schema) || alignmentMetric.matrix.same(ProteinSubstitutionMatrix.BLOSUM62),
				"Custom scores require a matrix-bound schema; legacy headers cannot discard score provenance");
			this.spectrumWidth=monomerSchema(schema) ? MONOMER_WIDTH : SPECTRUM_WIDTH;
			require(k>=4 && k<=6 && validDdlWidth(width),"Reference requires K4..K6 and 512 or 1024 registers");
			require(totalFamilies>=0 && liveFamilies>=0 && emptyFamilies>=0 &&
				(long)liveFamilies+emptyFamilies==totalFamilies,"Reference header census mismatch");
			this.schema=schema; this.k=k; this.width=width; this.seed=seed; this.totalFamilies=totalFamilies;
			this.liveFamilies=liveFamilies; this.emptyFamilies=emptyFamilies; this.formulaVersion=formulaVersion;
			this.spectrumEncoding=spectrumEncoding;
		}
		public final String schema,formulaVersion,spectrumEncoding;
		public final int k,width,totalFamilies,liveFamilies,emptyFamilies,spectrumWidth;
		public final long seed;
		public final double scalarZ;
		public final FamilyAlignmentMetric alignmentMetric;
	}

	/** Legacy consumers must fail rather than silently apply historical arithmetic to new gates. */
	public static Header parseHeader(byte[] line,LineParser1 lp){
		final Header header=parseMetricBoundHeader(line,lp);
		require(header.alignmentMetric.same(FamilyAlignmentMetric.LEGACY),
			"This reference uses a bound alignment metric; the consumer must explicitly support metric-aware loading");
		return header;
	}
	/** Parses either format for consumers that actually honor Header.alignmentMetric. */
	public static Header parseMetricBoundHeader(byte[] line,LineParser1 lp){
		require(line!=null && lp!=null,"Missing header/parser"); lp.set(line);
		require((lp.terms()==11 || lp.terms()==13 || lp.terms()==14 || lp.terms()==15) &&
			lp.termEquals("#family_reference",0),"Expected legacy, metric-bound or matrix-bound reference header");
		final String schema=lp.parseString(1);
		final boolean bound=metricSchema(schema);
		require(lp.terms()==(pairedSchema(schema) ? 15 : matrixSchema(schema) ? 14 : bound ? 13 : 11),
			"Alignment and matrix metadata must match the reference schema");
		FamilyAlignmentMetric metric=FamilyAlignmentMetric.LEGACY;
		if(bound){
			final FamilyAlignmentMetric.Mode mode=metricMode(lp.parseString(11));
			metric=new FamilyAlignmentMetric(mode,integer(lp,12),matrixSchema(schema) ?
				ProteinSubstitutionMatrix.parseBinding(lp.parseString(13)) : ProteinSubstitutionMatrix.BLOSUM62,
				pairedSchema(schema) ? integer(lp,14) : 0);
		}
		return new Header(schema,integer(lp,2),integer(lp,3),signedLong(lp,4),integer(lp,5),integer(lp,6),
			integer(lp,7),lp.parseString(8),lp.parseString(9),finiteDouble(lp,10),metric);
	}
	/** The header occupies one full five-line record, keeping batch boundaries aligned. */
	public static Header parseHeader(byte[][] lines,int offset,LineParser1 lp){
		return parseHeader(lines,offset,lp,false);
	}
	/** Five-line counterpart for metric-aware consumers; payload layout remains unchanged. */
	public static Header parseMetricBoundHeader(byte[][] lines,int offset,LineParser1 lp){
		return parseHeader(lines,offset,lp,true);
	}
	private static Header parseHeader(byte[][] lines,int offset,LineParser1 lp,boolean metricAware){
		require(lines!=null && offset>=0 && offset<=lines.length-LINES_PER_FAMILY,"Incomplete five-line reference header");
		final Header header=metricAware ? parseMetricBoundHeader(lines[offset],lp) : parseHeader(lines[offset],lp);
		for(int i=0;i<HEADER_LABELS.length;i++){
			final byte[] label=i==0 && alignmentStateSchema(header.schema) ? ALIGNMENT_GATES_LABEL :
				i==2 && header.width==512 ? COMPACT_DDL_LABEL :
				i==3 && header.spectrumWidth==MONOMER_WIDTH ? MONOMER_LABEL : HEADER_LABELS[i];
			require(Arrays.equals(lines[offset+i+1],label),"Invalid reference header label line "+(i+2));
		}
		return header;
	}
	public static ByteBuilder appendHeader(Header h,ByteBuilder out){
		require(h!=null && out!=null,"Missing header/output");
		out.append("#family_reference").tab().append(h.schema).tab().append(h.k).tab().append(h.width).tab()
			.append(h.seed).tab()
			.append(h.totalFamilies).tab().append(h.liveFamilies).tab().append(h.emptyFamilies).tab()
			.append(h.formulaVersion).tab().append(h.spectrumEncoding)
			.tab().append(Double.toString(h.scalarZ));
		if(metricSchema(h.schema)){out.tab().append(h.alignmentMetric.name()).tab().append(h.alignmentMetric.minimumSpan);}
		if(matrixSchema(h.schema)){out.tab().append(h.alignmentMetric.matrix.binding());}
		if(pairedSchema(h.schema)){out.tab().append(h.alignmentMetric.minimumPaired);}
		out.nl();
		for(int i=0;i<HEADER_LABELS.length;i++){
			out.append(i==0 && alignmentStateSchema(h.schema) ? ALIGNMENT_GATES_LABEL :
				i==2 && h.width==512 ? COMPACT_DDL_LABEL :
				i==3 && h.spectrumWidth==MONOMER_WIDTH ? MONOMER_LABEL : HEADER_LABELS[i]).nl();
		}
		return out;
	}
	/** Reference widths are explicit on the wire; the legacy default remains 1024. */
	static boolean validDdlWidth(final int width){return width==512 || width==DDL_WIDTH;}
	static boolean alignmentStateSchema(String schema){
		return ALIGNMENT_SCHEMA.equals(schema) || ALIGNMENT_MONOMER_SCHEMA.equals(schema) || matrixSchema(schema);
	}
	static boolean matrixSchema(String schema){
		return MATRIX_SCHEMA.equals(schema) || MATRIX_MONOMER_SCHEMA.equals(schema) || pairedSchema(schema);
	}
	private static boolean pairedSchema(String schema){
		return PAIRED_SCHEMA.equals(schema) || PAIRED_MONOMER_SCHEMA.equals(schema);
	}
	private static boolean pinnedMetric(FamilyAlignmentMetric.Mode mode){
		return mode==FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY ||
			mode==FamilyAlignmentMetric.Mode.PIN_LEFT_TRIM_RIGHT_ALIGNED_QUERY ||
			mode==FamilyAlignmentMetric.Mode.PIN_RIGHT_OVERHANG_ONLY_ALIGNED_QUERY;
	}
	private static FamilyAlignmentMetric.Mode metricMode(String label){
		for(FamilyAlignmentMetric.Mode mode:FamilyAlignmentMetric.Mode.values()){
			if(mode.label.equals(label)){
				require(pinnedMetric(mode),"Metric-bound references require a pinned aligned-query contract");
				return mode;
			}
		}
		throw new IllegalArgumentException("Unknown family alignment metric: "+label);
	}
	private static boolean metricSchema(String schema){
		return METRIC_SCHEMA.equals(schema) || METRIC_MONOMER_SCHEMA.equals(schema) || alignmentStateSchema(schema);
	}
	private static boolean monomerSchema(String schema){
		return MONOMER_SCHEMA.equals(schema) || METRIC_MONOMER_SCHEMA.equals(schema) ||
			ALIGNMENT_MONOMER_SCHEMA.equals(schema) || MATRIX_MONOMER_SCHEMA.equals(schema) ||
			PAIRED_MONOMER_SCHEMA.equals(schema);
	}

	public static FamilyReference parseFiveLines(byte[][] lines,int offset,LineParser1 lp){
		return parseFiveLines(lines,offset,lp,null);
	}

	/** Optional caller-owned worker timing; null performs no clock reads.
	 * Only successful records contribute. Intervals are disjoint elapsed time, not CPU time.
	 * DTO construction includes its validation and norm calculation, not the norm alone.
	 */
	public static FamilyReference parseFiveLines(byte[][] lines,int offset,LineParser1 lp,ParseTiming timing){
		return parseFiveLines(lines,offset,lp,timing,false,false,false);
	}

	/** Worker batches may precede the header. Final census must enforce both DDL and spectrum widths. */
	public static FamilyReference parseFiveLinesAnySpectrum(byte[][] lines,int offset,LineParser1 lp,ParseTiming timing){
		return parseFiveLines(lines,offset,lp,timing,true,true,true);
	}

	private static FamilyReference parseFiveLines(byte[][] lines,int offset,LineParser1 lp,ParseTiming timing,
			boolean allowMonomers,boolean allowAlignmentState,boolean allowCompactDdl){
		final long start=timing==null ? 0 : System.nanoTime();
		require(lines!=null && lp!=null && offset>=0 && offset<=lines.length-LINES_PER_FAMILY,"Incomplete five-line record");
		for(int i=0;i<LINES_PER_FAMILY;i++){require(lines[offset+i]!=null,"Null family payload line");}
		lp.set(lines[offset]); require(lp.terms()==6 && lp.termEquals("#family",0),"Expected six-field #family metadata");
		final String id=lp.parseString(1); final int rank=integer(lp,2);
		final long members=nonnegativeLong(lp,3),residues=nonnegativeLong(lp,4),words=nonnegativeLong(lp,5);
		if(members==0){
			for(int i=1;i<5;i++){
				final byte[] x=lines[offset+i];
				require(x.length==2 && x[0]=='N' && x[1]=='A',
					"Empty family requires four literal NA payload lines");
			}
			final long metadataEnd=timing==null ? 0 : System.nanoTime();
			require(residues==0 && words==0,"Empty family must have zero residue/word counts");
			final FamilyReference result=FamilyReference.empty(id,rank);
			if(timing!=null){
				final long end=System.nanoTime();
				timing.records++; timing.metadataGatesNanos+=metadataEnd-start; timing.dtoNormNanos+=end-metadataEnd;
			}
			return result;
		}
		lp.set(lines[offset+1]);
		final int gateTerms=lp.terms();
		final boolean explicitAlignment=gateTerms==18;
		require(gateTerms==17 || (allowAlignmentState && explicitAlignment),
			"Expected six A48 gates, scalar mask, optional alignment state and ten decimal bounds");
		final float cl=FamilyGateCodec.readLower(lp,0),il=FamilyGateCodec.readLower(lp,1),bl=FamilyGateCodec.readLower(lp,2),
			cu=FamilyGateCodec.readCosineAccept(lp,3),iu=FamilyGateCodec.readUpper(lp,4),
			bf=FamilyGateCodec.readFinalBlosum(lp,5);
		final int scalarMask=integer(lp,6), scalarOffset=explicitAlignment ? 8 : 7;
		final FamilyReference.AlignmentCalibration alignmentCalibration=explicitAlignment ? alignmentCalibration(lp,7) :
			FamilyReference.AlignmentCalibration.ENABLED;
		final double[] scalarBounds=new double[10];
		for(int i=0;i<scalarBounds.length;i++){scalarBounds[i]=finiteDouble(lp,scalarOffset+i);}
		final long metadataEnd=timing==null ? 0 : System.nanoTime();
		final byte[] raw=lines[offset+2]; final ProteinSequence consensus=new ProteinSequence(id,raw);
		lp.set(lines[offset+3]); final int ddlWidth=lp.terms();
		require(ddlWidth==DDL_WIDTH || allowCompactDdl && ddlWidth==512,
			"Expected 1024 A48 DDL registers or explicitly supported 512-register payload");
		final char[] ddl=new char[ddlWidth];
		for(int i=0;i<ddl.length;i++){ddl[i]=(char)a48(lp,i,65535L);}
		final long consensusDdlEnd=timing==null ? 0 : System.nanoTime();
		lp.set(lines[offset+4]); final int spectrumWidth=lp.terms();
		require(spectrumWidth==SPECTRUM_WIDTH || allowMonomers && spectrumWidth==MONOMER_WIDTH,
			"Expected 400 A48 dimer values or explicitly supported 20 projected monomers");
		final float[] spectrum=new float[spectrumWidth];
		for(int i=0;i<spectrum.length;i++){spectrum[i]=Float.intBitsToFloat((int)a48(lp,i,0xffffffffL));}
		final long spectrumEnd=timing==null ? 0 : System.nanoTime();
		final FamilyReference result=explicitAlignment ?
			new FamilyReference(id,rank,members,residues,words,consensus,ddl,spectrum,cl,il,bl,cu,iu,bf,raw,
				scalarMask,scalarBounds,alignmentCalibration) :
			new FamilyReference(id,rank,members,residues,words,consensus,ddl,spectrum,cl,il,bl,cu,iu,bf,raw,
				scalarMask,scalarBounds);
		if(timing!=null){
			final long end=System.nanoTime();
			timing.records++; timing.metadataGatesNanos+=metadataEnd-start;
			timing.consensusDdlNanos+=consensusDdlEnd-metadataEnd;
			timing.spectrumNanos+=spectrumEnd-consensusDdlEnd; timing.dtoNormNanos+=end-spectrumEnd;
		}
		return result;
	}

	/** Mutable worker-local sums; merge only after the source worker has joined.
	 * Not thread-safe. Empty records contribute metadata/NA validation and DTO work only.
	 * Summed worker elapsed intervals overlap across workers and are not loader wall time.
	 */
	public static final class ParseTiming{
		public void add(ParseTiming other){
			require(other!=null,"Missing completed worker timing to merge");
			records+=other.records; metadataGatesNanos+=other.metadataGatesNanos;
			consensusDdlNanos+=other.consensusDdlNanos; spectrumNanos+=other.spectrumNanos;
			dtoNormNanos+=other.dtoNormNanos;
		}
		public long records,metadataGatesNanos,consensusDdlNanos,spectrumNanos,dtoNormNanos;
	}

	public static ByteBuilder appendFiveLines(FamilyReference f,ByteBuilder out){
		require(f!=null && (f.members==0 || f.ddl.length==DDL_WIDTH),
			"512-register output requires an explicit matching reference header");
		return appendFiveLines(f,out,false);
	}

	/** New stateful schemas must opt in so unavailable alignment gates cannot be flattened to legacy bytes. */
	public static ByteBuilder appendFiveLines(FamilyReference f,ByteBuilder out,Header header){
		require(header!=null && f!=null,"Missing family/reference header for payload output");
		require(f.members==0 || f.ddl.length==header.width && f.spectrum.length==header.spectrumWidth,
			"Reference payload widths disagree with output header: "+f.id);
		return appendFiveLines(f,out,alignmentStateSchema(header.schema));
	}

	private static ByteBuilder appendFiveLines(FamilyReference f,ByteBuilder out,boolean explicitAlignment){
		require(f!=null && out!=null,"Missing family/output");
		require(f.members==0 || f.explicitAlignmentCalibration==explicitAlignment,
			"Reference alignment-state row disagrees with package header: "+f.id);
		out.append("#family").tab().append(f.id).tab().append(f.rank).tab().append(f.members).tab()
			.append(f.residues).tab().append(f.validWords).nl();
		if(f.members==0){return out.append("NA\nNA\nNA\nNA\n");}
		FamilyGateCodec.appendLower(out,f.cosineLower).tab(); FamilyGateCodec.appendLower(out,f.identityLower).tab();
		FamilyGateCodec.appendLower(out,f.blosumLower).tab(); FamilyGateCodec.appendCosineAccept(out,f.cosineAccept).tab();
		FamilyGateCodec.appendUpper(out,f.identityUpper).tab();
		FamilyGateCodec.appendFinalBlosum(out,f.blosumFinal).tab().append(f.scalarMask);
		if(explicitAlignment){out.tab().append(f.alignmentCalibration.label);}
		for(double bound:f.scalarBounds){out.tab().append(Double.toString(bound));} out.nl();
		out.append(f.consensusRaw).nl();
		for(int i=0;i<f.ddl.length;i++){if(i>0){out.tab();} out.appendA48((long)f.ddl[i]);} out.nl();
		for(int i=0;i<f.spectrum.length;i++){if(i>0){out.tab();} out.appendFloatA48(f.spectrum[i]);}
		return out.nl();
	}

	private static double finiteDouble(LineParser1 lp,int term){
		// Preserve exact decimal-double roundtrips, including exponent notation;
		// Parse's fast decimal path does not support exponents.
		final double value=Double.parseDouble(lp.parseString(term));
		require(Double.isFinite(value),"Scalar fields must be finite decimal doubles"); return value;
	}
	private static long a48(LineParser1 lp,int term,long max){
		final int length=lp.setBounds(term),start=lp.a(); final byte[] line=lp.line();
		require(length>0 && length<=6 && (length==1 || line[start]!='0'),"Noncanonical A48 token");
		for(int i=start;i<lp.b();i++){require(line[i]>=48 && line[i]<=111,"Invalid A48 alphabet");}
		final long value=lp.parseLongA48(term); require(value>=0 && value<=max,"A48 value exceeds field range"); return value;
	}
	private static int integer(LineParser1 lp,int term){
		final long x=nonnegativeLong(lp,term);
		require(x<=Integer.MAX_VALUE,"Reference integer overflow");
		return (int)x;
	}
	private static FamilyReference.AlignmentCalibration alignmentCalibration(LineParser1 lp,int term){
		if(lp.termEquals(FamilyReference.AlignmentCalibration.ENABLED.label,term)){
			return FamilyReference.AlignmentCalibration.ENABLED;
		}
		if(lp.termEquals(FamilyReference.AlignmentCalibration.UNAVAILABLE.label,term)){
			return FamilyReference.AlignmentCalibration.UNAVAILABLE;
		}
		throw new IllegalArgumentException("Unknown alignment calibration availability");
	}
	private static long nonnegativeLong(LineParser1 lp,int term){
		final long x=signedLong(lp,term);
		require(x>=0,"Negative reference count/rank");
		return x;
	}
	private static long signedLong(LineParser1 lp,int term){
		final int length=lp.setBounds(term); require(length>0,"Empty decimal integer");
		final byte[] bytes=lp.line(); int pos=lp.a(); final boolean negative=bytes[pos]=='-'; if(negative){pos++;}
		require(pos<lp.b(),"Missing decimal digits");
		long value=0; final long limit=negative ? Long.MIN_VALUE : -Long.MAX_VALUE;
		for(;pos<lp.b();pos++){
			final int digit=bytes[pos]-'0';
			require(digit>=0 && digit<=9 && value>=limit/10,"Invalid or overflowing decimal integer");
			value*=10; require(value>=limit+digit,"Decimal integer overflow"); value-=digit;
		}
		return negative ? value : -value;
	}
	private static void require(boolean condition,String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}
}
