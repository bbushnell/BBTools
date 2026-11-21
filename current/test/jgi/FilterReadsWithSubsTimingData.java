package test.jgi;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.Parse;
import structures.ByteBuilder;

/** Generates and validates the nonterminal eight-template FilterReadsWithSubs timing SAM.
 * This is data preparation only; it does not run FilterReadsWithSubs or time old/fixed code.
 * @author Jean, Shinobu
 */
public final class FilterReadsWithSubsTimingData{

	/** Runs mode=generate or mode=check for exactly one cycle count. */
	public static void main(final String[] args) throws Exception{
		String mode=null, in=null, out=null, outdir=null;
		int cycles=-1;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("mode")){mode=pair[1];}
			else if(pair[0].equals("in")){in=pair[1];}
			else if(pair[0].equals("out")){out=pair[1];}
			else if(pair[0].equals("outdir")){outdir=pair[1];}
			else if(pair[0].equals("cycles")){cycles=Parse.parseIntKMG(pair[1]);}
			else{throw new IllegalArgumentException("Unknown timing-data option: "+arg);}
		}
		require("generate".equals(mode) || "check".equals(mode),
			"Set mode=generate or mode=check.");
		require(cycles>0, "Set cycles to a positive integer; cycles="+cycles);
		require(outdir!=null && !outdir.isEmpty(), "Set a fresh outdir.");
		final Path root=Paths.get(outdir).toAbsolutePath().normalize();
		Files.createDirectory(root);
		final FilterReadsWithSubsTimingData helper=
			new FilterReadsWithSubsTimingData(root, cycles);
		if("generate".equals(mode)){
			require(in==null, "mode=generate writes out=<sam> and takes no in=<sam>.");
			require(out!=null && !out.isEmpty(), "mode=generate requires out=<sam>.");
			helper.generate(out);
		}else{
			require(out==null, "mode=check reads in=<sam> and takes no out=<sam>.");
			require(in!=null && !in.isEmpty(), "mode=check requires in=<sam>.");
			helper.check(in, "retained");
		}
	}

	/** Prepares the output directory and template invariants. */
	private FilterReadsWithSubsTimingData(final Path root_, final int cycles_){
		root=root_;
		cycles=cycles_;
		expectedRecords=Math.multiplyExact(8, cycles);
		expectedBases=Math.multiplyExact(800L, cycles);
		assert(TEMPLATES.length==8) : "The JGI154 timing proposal in plans/jgi-maintenance.md has exactly eight templates.";
		for(final Template template : TEMPLATES){
			assert(template.queryBases==100) : "Each JGI154 timing template must have 100 query bases: "+
				template.name+" has "+template.queryBases;
			assert(template.referenceSpan<=102) : "With POS<=1000 and LN=2000, every agreed template endpoint must be <=1101; "+
				template.name+" span="+template.referenceSpan;
		}
	}

	/** Writes the exact SAM input and then independently validates it. */
	private void generate(final String fname) throws IOException{
		final File f=new File(fname);
		require(!f.exists(), "Refusing to overwrite generated SAM: "+fname);
		ByteStreamWriter writer=null;
		boolean error=false;
		try{
			writer=new ByteStreamWriter(fname, false, false, false);
			writer.start();
			final ByteBuilder bb=new ByteBuilder(1<<20);
			bb.append("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:2000\n");
			for(int i=0; i<expectedRecords; i++){
				appendRecord(bb, i);
				if(bb.length()>(1<<19)){writer.print(bb); bb.clear();}
			}
			if(!bb.isEmpty()){writer.print(bb);}
		}finally{
			if(writer!=null){error=writer.poisonAndWait();}
		}
		require(!error, "Generated SAM writer reported an error for "+fname);
		check(fname, "generated");
	}

	/** Validates the generated or retained SAM records against ordinal-derived expectations. */
	private void check(final String fname, final String label) throws IOException{
		final File f=new File(fname);
		require(f.isFile(), "Expected readable SAM file: "+fname);
		final Result result=new Result();
		final int[] starts=new int[MAX_FIELDS], stops=new int[MAX_FIELDS];
		ByteFile bf=null;
		boolean error=false;
		try{
			bf=ByteFile.makeByteFile(fname, false);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length<1){throw new IllegalArgumentException("Blank SAM line in "+fname);}
				if(line[0]=='@'){
					validateHeader(line, result);
				}else{
					if(result.records>=expectedRecords){
						throw new IllegalArgumentException("SAM contains more than "+expectedRecords+" records.");
					}
					validateRecord(line, result.records, starts, stops);
					result.records++;
					result.bases+=READ_LENGTH;
				}
			}
		}finally{
			if(bf!=null){error=bf.close();}
		}
		require(!error, "SAM reader reported an error for "+fname);
		require(result.headers==2, "Expected 2 ordered headers but observed "+result.headers);
		require(result.records==expectedRecords, "Expected "+expectedRecords+
			" retained records from "+cycles+" cycles but saw "+result.records);
		require(result.bases==expectedBases, "Expected "+expectedBases+
			" retained query bases from "+cycles+" cycles but saw "+result.bases);
		writeSummary(label, result);
	}

	/** Appends the deterministic record with ordinal-specific QNAME and POS. */
	private static void appendRecord(final ByteBuilder bb, final int ordinal){
		final Template template=TEMPLATES[ordinal&7];
		bb.append("bench_").append(ordinal).tab();
		bb.append('0').tab();
		bb.append("ref\t");
		bb.append(1+(ordinal%1000)).tab();
		bb.append("60\t");
		bb.append(template.cigar).tab();
		bb.append("*\t0\t0\t");
		bb.append(C_BASES).tab();
		bb.append(I_QUALS).tab();
		bb.append("MD:Z:").append(template.md).nl();
	}

	/** Checks one expected record without trusting any field from the input line. */
	private void validateRecord(final byte[] line, final int ordinal,
			final int[] starts, final int[] stops){
		final int fields=split(line, starts, stops);
		if(fields!=12){
			throw new IllegalArgumentException("SAM record "+ordinal+" must contain exactly 12 fields but had "+fields);
		}
		final Template template=TEMPLATES[ordinal&7];
		requireBench(line, starts[0], stops[0], ordinal);
		requireField(line, starts[1], stops[1], "0", ordinal, "FLAG");
		requireField(line, starts[2], stops[2], "ref", ordinal, "RNAME");
		requireInt(line, starts[3], stops[3], 1+(ordinal%1000), ordinal, "POS");
		requireField(line, starts[4], stops[4], "60", ordinal, "MAPQ");
		requireField(line, starts[5], stops[5], template.cigar, ordinal, "CIGAR");
		requireField(line, starts[6], stops[6], "*", ordinal, "RNEXT");
		requireField(line, starts[7], stops[7], "0", ordinal, "PNEXT");
		requireField(line, starts[8], stops[8], "0", ordinal, "TLEN");
		requireBytes(line, starts[9], stops[9], C_BASES, ordinal, "SEQ");
		requireBytes(line, starts[10], stops[10], I_QUALS, ordinal, "QUAL");
		if(!startsWith(line, starts[11], stops[11], "MD:Z:")){
			throw new IllegalArgumentException("SAM record "+ordinal+" is missing the MD:Z: prefix.");
		}
		requireField(line, starts[11]+5, stops[11], template.md, ordinal, "MD");
	}

	/** Validates the exact two generated or retained SAM headers. */
	private static void validateHeader(final byte[] line, final Result result){
		require(result.records==0, "SAM header appeared after records.");
		if(result.headers==0){
			requireField(line, 0, line.length, "@HD\tVN:1.6\tSO:unsorted", -1, "@HD");
		}else if(result.headers==1){
			requireField(line, 0, line.length, "@SQ\tSN:ref\tLN:2000", -1, "@SQ");
		}else{
			throw new IllegalArgumentException("Unexpected extra SAM header line.");
		}
		result.headers++;
	}

	/** Writes a compact validation receipt for root's timed launch recipes. */
	private void writeSummary(final String label, final Result result) throws IOException{
		final Path path=root.resolve(label+"-summary.txt");
		ByteStreamWriter writer=null;
		boolean error=false;
		try{
			writer=new ByteStreamWriter(path.toString(), false, false, false);
			writer.start();
			final ByteBuilder bb=new ByteBuilder(256);
			bb.append("FILTERSUBS_TIMING_DATA_PASS label=").append(label);
			bb.append(" cycles=").append(cycles);
			bb.append(" records=").append(result.records);
			bb.append(" expectedRecords=").append(expectedRecords);
			bb.append(" bases=").append(result.bases);
			bb.append(" expectedBases=").append(expectedBases);
			bb.append(" deletionTemplates=4 deletionOpsPerCycle=5 failures=0\n");
			writer.print(bb);
			System.out.print(bb.toString());
		}finally{
			if(writer!=null){error=writer.poisonAndWait();}
		}
		require(!error, "Summary writer reported an error for "+path);
	}

	/** Splits a tab-delimited SAM row into caller-supplied start/stop arrays. */
	private static int split(final byte[] line, final int[] starts, final int[] stops){
		int fields=0, start=0;
		for(int i=0; i<=line.length; i++){
			if(i==line.length || line[i]=='\t'){
				if(fields>=starts.length){
					throw new IllegalArgumentException("SAM record has too many fields; max="+starts.length);
				}
				starts[fields]=start;
				stops[fields]=i;
				fields++;
				start=i+1;
			}
		}
		return fields;
	}

	/** Checks bench_i without allocating a String per record. */
	private static void requireBench(final byte[] line, final int start, final int stop,
			final int ordinal){
		if(!startsWith(line, start, stop, "bench_")){
			throw new IllegalArgumentException("SAM record "+ordinal+
				" has an unexpected QNAME prefix.");
		}
		requireDecimal(line, start+6, stop, ordinal, ordinal, "QNAME");
	}

	/** Checks one canonical decimal field against an expected value. */
	private static void requireInt(final byte[] line, final int start, final int stop,
			final int expected, final int ordinal, final String name){
		requireDecimal(line, start, stop, expected, ordinal, name);
	}

	/** Checks one ASCII field against an expected literal. */
	private static void requireField(final byte[] line, final int start, final int stop,
			final String expected, final int ordinal, final String name){
		if(!equals(line, start, stop, expected)){
			throw new IllegalArgumentException("SAM record "+ordinal+
				" has unexpected "+name+"; expected "+expected);
		}
	}

	/** Checks one byte field against an expected repeated payload. */
	private static void requireBytes(final byte[] line, final int start, final int stop,
			final byte[] expected, final int ordinal, final String name){
		if(stop-start!=expected.length){
			throw new IllegalArgumentException("SAM record "+ordinal+" "+name+
				" length expected "+expected.length+" but observed "+(stop-start));
		}
		for(int i=0; i<expected.length; i++){
			if(line[start+i]!=expected[i]){
				throw new IllegalArgumentException("SAM record "+ordinal+" "+name+
					" differs at offset "+i);
			}
		}
	}

	/** Checks a field's bytes against Java's canonical nonnegative decimal spelling. */
	private static void requireDecimal(final byte[] line, final int start, final int stop,
			final int expected, final int ordinal, final String name){
		assert(expected>=0) : "This fixture emits only nonnegative decimal fields; "+name+"="+expected;
		int digits=(expected==0 ? 1 : 0);
		for(int x=expected; x>0; x/=10){digits++;}
		if(stop-start!=digits){throw new IllegalArgumentException("SAM record "+ordinal+
			" "+name+" is not the canonical decimal "+expected);}
		int divisor=1;
		for(int i=1; i<digits; i++){divisor*=10;}
		int x=expected;
		for(int i=start; i<stop; i++){
			final byte expectedByte=(byte)('0'+(x/divisor));
			if(line[i]!=expectedByte){throw new IllegalArgumentException("SAM record "+ordinal+
				" "+name+" is not the canonical decimal "+expected);}
			x%=divisor;
			divisor/=10;
		}
	}

	/** Tests whether a field starts with a literal. */
	private static boolean startsWith(final byte[] line, final int start, final int stop,
			final String expected){
		final int len=expected.length();
		if(stop-start<len){return false;}
		for(int i=0; i<len; i++){
			if(line[start+i]!=(byte)expected.charAt(i)){return false;}
		}
		return true;
	}

	/** Tests one whole field against a literal. */
	private static boolean equals(final byte[] line, final int start, final int stop,
			final String expected){
		if(stop-start!=expected.length()){return false;}
		return startsWith(line, start, stop, expected);
	}

	/** Repeats one ASCII byte. */
	private static byte[] repeat(final byte value, final int count){
		final byte[] array=new byte[count];
		for(int i=0; i<count; i++){array[i]=value;}
		return array;
	}

	/** Fails loudly for fixture, input or output mismatches. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** One exact CIGAR/MD row from the final JGI154 timing proposal. */
	private static final class Template{
		Template(final String name_, final String cigar_, final String md_,
				final int queryBases_, final int referenceSpan_){
			name=name_;
			cigar=cigar_;
			md=md_;
			queryBases=queryBases_;
			referenceSpan=referenceSpan_;
		}

		final String name;
		final String cigar;
		final String md;
		final int queryBases;
		final int referenceSpan;
	}

	/** Accumulates one SAM validation pass. */
	private static final class Result{
		int headers;
		int records;
		long bases;
	}

	private final Path root;
	private final int cycles;
	private final int expectedRecords;
	private final long expectedBases;

	private static final int READ_LENGTH=100;
	private static final int MAX_FIELDS=64;
	private static final byte[] C_BASES=repeat((byte)'C', READ_LENGTH);
	private static final byte[] I_QUALS=repeat((byte)'I', READ_LENGTH);
	private static final Template[] TEMPLATES={
		new Template("perfect", "100=", "100", 100, 100),
		new Template("substitution", "100M", "50A49", 100, 100),
		new Template("internal_deletion", "40=1D60=", "40^A60", 100, 101),
		new Template("internal_insertion", "40=1I59=", "99", 100, 99),
		new Template("two_internal_deletions", "25=1D25=1D50=", "25^A25^C50", 100, 102),
		new Template("substitution_insertion", "30M1I69M", "10A88", 100, 99),
		new Template("leading_clip_deletion", "5S45=1D50=", "45^G50", 100, 96),
		new Template("trailing_clip_deletion", "45=1D50=5S", "45^T50", 100, 96)
	};
}
