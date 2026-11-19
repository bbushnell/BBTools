package prot.cluster;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntList;

/** Joins independent structural probes to their matched unchanged reassignment control.
 * This reducer reports counts, not fitted coefficients or biological acceptance.
 * Raw proposal features remain in the producer table, joined by candidate_index.
 * Input provenance and identical execution settings remain the controller's duty.
 * @author Brian Bushnell, Fischl
 */
public final class DdlClustProbeReport{

	public static void main(final String[] args) throws Exception{
		final DdlClustProbeReport tool=new DdlClustProbeReport(args);
		try{tool.process();}finally{Shared.closeStream(tool.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses native config/argument handling; the count pattern has one candidate-ID placeholder. */
	public DdlClustProbeReport(String[] args){
		final PreParser pre=new PreParser(args,getClass(),false);
		args=pre.args; outstream=pre.outstream;
		final Parser parser=new Parser(); parser.overwrite=false;
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			final String key=pair[0].toLowerCase(Locale.ROOT);
			final String value=pair.length>1 && !pair[1].equalsIgnoreCase("null") ? pair[1] : null;
			if(key.equals("results")){results=value;}
			else if(key.equals("control")){control=value;}
			else if(key.equals("counts")){pattern=value;}
			else if(key.equals("expectedinputs")){expectedInputs=Parse.parseKMG(value);}
			else if(!parser.parse(arg,key,value)){throw new IllegalArgumentException("Unknown probe-report parameter: "+arg);}
		}
		out=parser.out1; overwrite=parser.overwrite;
		require(!parser.append,"Probe reports cannot append");
		require(results!=null && control!=null && pattern!=null && out!=null && expectedInputs>0,
			"Require results, control, counts, out and positive expectedinputs");
		require(pattern.indexOf('#')>=0 && pattern.indexOf('#')==pattern.lastIndexOf('#'),
			"Counts pattern must contain exactly one # for candidate_index");
	}

	/** Validates every source and computes every response before opening the report. */
	public void process() throws IOException{
		final ArrayList<Path> paths=new ArrayList<Path>();
		register(paths,results,true); register(paths,control,true);
		final ArrayList<Probe> probes=readProbes();
		final IntList originalIds=new IntList();
		final IntHashMap originals=new IntHashMap();
		for(Probe p:probes){
			addId(originalIds,originals,p.a);
			if(p.b>=0){addId(originalIds,originals,p.b);}
			p.countPath=pattern.replace("#",Integer.toString(p.id));
			if(!p.feasible){
				require(!Files.exists(Paths.get(p.countPath),LinkOption.NOFOLLOW_LINKS),
					"Rejected candidate has a stale count path: "+p.countPath);
			}
			register(paths,p.countPath,p.feasible);
		}
		register(paths,out,false);
		require(overwrite || !Files.exists(Paths.get(out),LinkOption.NOFOLLOW_LINKS),"Report exists and overwrite=false: "+out);
		require(Tools.testOutputFiles(overwrite,false,false,out),"Probe-report output is unavailable");
		final Counts baseline=readCounts(control,Arrays.copyOf(originalIds.array,originalIds.size),null);
		for(Probe p:probes){
			p.control=baseline.count(p.a);
			if(p.b>=0){p.control=Math.addExact(p.control,baseline.count(p.b));}
			if(p.feasible){
				final Counts modified=readCounts(p.countPath,p.replacements,baseline);
				p.modified=modified.sum;
			}
		}
		write(probes);
		outstream.println("DDLCLUST_PROBE_REPORT_PASS candidates="+probes.size()+" input_records="+expectedInputs);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Input accounting      ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads retained proposal rows in producer order, including genuine counts on early rejection. */
	private ArrayList<Probe> readProbes(){
		final ArrayList<Probe> probes=new ArrayList<Probe>();
		final IntHashMap ids=new IntHashMap();
		final ByteFile reader=ByteFile.makeByteFile(results,false);
		try{
			header(reader.nextLine(),StructuralProbeBatch.RESULT_HEADER,"producer");
			final LineParser1 row=new LineParser1('\t');
			final ByteBuilder removals=new ByteBuilder();
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				row.set(line); require(row.terms()==31,"Producer rows require exactly 31 columns");
				final int id=identifier(row,0,false),a=identifier(row,5,false),b=identifier(row,6,true);
				decimal(row,1); // Validate seed metadata, including its signed-long range.
				require(ids.get(id)<0,"Duplicate candidate_index"); ids.put(id,1);
				final boolean split=row.termEquals("fission",2);
				require(split || row.termEquals("fusion",2),"Unknown probe operation");
				final long feasible=nonnegative(row,3);
				require(feasible<=1,"Feasible must be 0 or 1");
				require(split ? b==-1 : b>a,"Probe endpoints disagree with native ordered operation");
				require(row.setBounds(4)>0 && !row.termEquals("NA",4),"Missing producer reason");
				removals.clear().append(a); if(b>=0){removals.append(',').append(b);}
				require(equalTerm(row,7,removals),"Removal IDs disagree with original endpoints");
				require(row.setBounds(8)>0 && !row.termEquals("NA",8),"Missing overlay path");
				final long initial=nonnegative(row,30);
				require(initial>0 && initial<=expectedInputs,"Initial affected members exceed the input census or are empty");
				nonnegative(row,17); finite(row,18);
				final int[] replacements;
				if(feasible==1){
					final int parent=identifier(row,9,true),left=identifier(row,10,false),
						right=identifier(row,11,false),union=identifier(row,12,true);
					final long parentN=nonnegative(row,13),leftN=nonnegative(row,14),
						rightN=nonnegative(row,15),unionN=nonnegative(row,16);
					require(left!=right && leftN>0 && rightN>0 && Math.addExact(leftN,rightN)==initial,
						"Proposal child/endpoint sizes must conserve initial affected members");
					if(split){
						require(parent==a && union==-1 && parentN==initial && unionN==0 && left==a && right!=a,
							"Fission IDs/counts disagree with its one-parent, two-child shape");
						replacements=new int[]{left,right};
					}else{
						require(parent==-1 && parentN==0 && union==a && unionN==initial && left==a && right==b,
							"Fusion IDs/counts disagree with its two-parent, one-union shape");
						replacements=new int[]{union};
					}
					for(int i=19; i<30; i++){finite(row,i);}
				}else{
					for(int i=9; i<17; i++){require(row.termEquals("NA",i),"Rejected proposal must not retain model IDs/counts");}
					for(int i=19; i<30; i++){require(row.termEquals("NA",i),"Rejected proposal must not retain feature values");}
					replacements=new int[0];
				}
				probes.add(new Probe(id,a,b,split,feasible==1,row.parseByteArray(4),initial,replacements));
			}
		}finally{close(reader,results);}
		require(!probes.isEmpty(),"Empty producer table");
		return probes;
	}

	/** Exact requested-ID coverage distinguishes an observed zero from an absent measurement.
	 * Each file has its own namespace: fresh model IDs may recur across independent probes.
	 */
	private Counts readCounts(final String path,final int[] requested,final Counts baseline){
		assert(requested.length>0) : "Every count table must describe at least one explicit affected family";
		final Counts counts=new Counts(requested);
		final ByteFile reader=ByteFile.makeByteFile(path,false);
		try{
			header(reader.nextLine(),COUNTS_HEADER,"affected counts");
			final LineParser1 row=new LineParser1('\t');
			final long[] totals=new long[5]; int seen=0;
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				row.set(line); require(row.terms()==7,"Affected-count rows require exactly seven columns");
				final int id=identifier(row,0,false),offset=counts.offsets.get(id);
				require(offset>=0 && counts.values[offset]<0,"Unexpected or duplicate affected family ID in "+path);
				final long value=nonnegative(row,1);
				for(int i=0; i<5; i++){totals[i]=nonnegative(row,i+2);}
				require(totals[0]==expectedInputs && totals[1]<=totals[0] && totals[2]<=totals[1] &&
					totals[3]==totals[0]-totals[2] && totals[4]==totals[0]-totals[1],"Inconsistent count-table census in "+path);
				if(seen==0){System.arraycopy(totals,0,counts.totals,0,5);}
				else{require(Arrays.equals(totals,counts.totals),"Count-table totals change between family rows in "+path);}
				require(value<=totals[2],"Affected family count exceeds total assigned in "+path);
				counts.values[offset]=value; counts.sum=Math.addExact(counts.sum,value); seen++;
			}
			require(seen==requested.length,"Missing affected family ID in "+path);
			require(counts.sum<=counts.totals[2],"Affected family sum exceeds total assigned in "+path);
			if(baseline!=null){
				require(counts.totals[0]==baseline.totals[0] && counts.totals[1]==baseline.totals[1] &&
					counts.totals[4]==baseline.totals[4],"Modified and control inputs/eligibility/filtering disagree");
			}
			return counts;
		}finally{close(reader,path);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	/** A rejected measurement is unavailable, while explicit zero remains a numeric observation. */
	private void write(final ArrayList<Probe> probes){
		final ByteStreamWriter writer=new ByteStreamWriter(out,overwrite,false,true); writer.start();
		try{
			writer.println(REPORT_HEADER);
			final ByteBuilder row=new ByteBuilder();
			for(Probe p:probes){
				row.clear().append(p.id).tab().append(p.split ? "fission" : "fusion").tab()
					.append(p.feasible ? 1 : 0).tab().append(p.reason).tab().append(p.a).tab().append(p.b)
					.tab().append(p.initial).tab().append(p.control).tab();
				if(p.feasible){row.append(p.modified).tab().append(Math.subtractExact(p.modified,p.initial));}
				else{row.append("NA\tNA");}
				row.tab().append(Math.subtractExact(p.control,p.initial)).tab();
				if(p.feasible){row.append(Math.subtractExact(p.modified,p.control));}else{row.append("NA");}
				row.tab().append(expectedInputs); writer.println(row);
			}
		}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Failed probe report output: "+out);}}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Parsing and paths      ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves existing symlinks and hardlinks; rejected candidate directories may not exist at all. */
	private static void register(final ArrayList<Path> paths,final String name,final boolean input) throws IOException{
		require(name!=null && !name.isEmpty() && name.indexOf('\t')<0 && name.indexOf('\n')<0 && name.indexOf('\r')<0,
			"Report paths must be nonempty unescaped TSV fields");
		final Path path=Paths.get(name).toAbsolutePath().normalize();
		require(path.getFileName()!=null && !FileFormat.isStdio(path.getFileName().toString()) &&
			!path.toString().equals("/dev/null"),"Report paths must name disk files");
		final Path real;
		if(input || Files.exists(path,LinkOption.NOFOLLOW_LINKS)){
			real=path.toRealPath(); require(Files.isRegularFile(real),"Report input/output must be a regular file");
		}else{
			Path ancestor=path.getParent();
			while(ancestor!=null && !Files.exists(ancestor,LinkOption.NOFOLLOW_LINKS)){ancestor=ancestor.getParent();}
			require(ancestor!=null,"Path has no existing ancestor");
			real=ancestor.toRealPath().resolve(ancestor.relativize(path)).normalize();
		}
		for(Path old:paths){
			require(!old.equals(real) && !(Files.exists(old) && Files.exists(real) && Files.isSameFile(old,real)),
				"Probe-report source/output paths alias");
		}
		paths.add(real);
	}
	/** Negative accumulation permits Long.MIN_VALUE seeds without permitting wrapped counts. */
	private static long decimal(final LineParser1 row,final int term){
		require(row.setBounds(term)>0,"Missing integer field");
		final byte[] line=row.line(); int at=row.a(); final int end=row.b();
		final boolean negative=line[at]=='-'; if(negative){at++;}
		require(at<end,"Missing integer digits");
		final long limit=negative ? Long.MIN_VALUE : -Long.MAX_VALUE; long value=0;
		for(; at<end; at++){
			final int digit=line[at]-'0';
			require(digit>=0 && digit<=9 && value>=limit/10,"Invalid or overflowing integer");
			value*=10; require(value>=limit+digit,"Overflowing integer"); value-=digit;
		}
		return negative ? value : -value;
	}
	private static long nonnegative(final LineParser1 row,final int term){
		final long value=decimal(row,term); require(value>=0,"Negative count field"); return value;
	}
	private static int identifier(final LineParser1 row,final int term,final boolean absent){
		final long value=decimal(row,term);
		require(value>=(absent ? -1 : 0) && value<Integer.MAX_VALUE,"Identifier out of range"); return (int)value;
	}
	private static void finite(final LineParser1 row,final int term){
		// Producer appendSlow emits exponents. This small metadata table needs the native exponent-capable parser.
		require(row.setBounds(term)>0 && Double.isFinite(Parse.parseDoubleSlow(row.line(),row.a(),row.b())),
			"Nonfinite producer numeric field");
	}
	private static boolean equalTerm(final LineParser1 parser,final int term,final ByteBuilder expected){
		if(parser.setBounds(term)!=expected.length()){return false;}
		for(int i=0; i<expected.length(); i++){if(parser.line()[parser.a()+i]!=expected.array[i]){return false;}}
		return true;
	}
	private static void addId(final IntList list,final IntHashMap offsets,final int id){
		if(offsets.get(id)<0){offsets.put(id,list.size); list.add(id);}
	}
	private static void header(final byte[] line,final String expected,final String description){
		require(line!=null && Arrays.equals(line,expected.getBytes(StandardCharsets.US_ASCII)),"Unexpected "+description+" header");
	}
	private static void close(final ByteFile reader,final String path){
		if(reader.close()){throw new IllegalStateException("Failed probe-report input: "+path);}
	}
	private static void require(final boolean valid,final String message){if(!valid){throw new IllegalArgumentException(message);}}

	/*--------------------------------------------------------------*/
	/*----------------           State              ----------------*/
	/*--------------------------------------------------------------*/

	/** Only compact result rows survive parsing; large sequence/model data are never loaded. */
	private static final class Probe{
		Probe(final int id_,final int a_,final int b_,final boolean split_,final boolean feasible_,
				final byte[] reason_,final long initial_,final int[] replacements_){
			id=id_; a=a_; b=b_; split=split_; feasible=feasible_; reason=reason_; initial=initial_; replacements=replacements_;
		}
		final int id,a,b;
		final boolean split,feasible;
		final byte[] reason;
		final long initial;
		final int[] replacements;
		String countPath;
		long control,modified;
	}
	/** Per-file ID lookup retains literal zero; -1 means the required row was not present. */
	private static final class Counts{
		Counts(final int[] ids){
			values=new long[ids.length]; Arrays.fill(values,-1);
			for(int i=0; i<ids.length; i++){
				require(offsets.get(ids[i])<0,"Repeated expected family ID"); offsets.put(ids[i],i);
			}
		}
		long count(final int id){
			final int offset=offsets.get(id);
			assert(offset>=0 && values[offset]>=0) : "Complete requested-ID coverage was checked before matched-control lookup";
			return values[offset];
		}
		final IntHashMap offsets=new IntHashMap();
		final long[] values,totals=new long[5];
		long sum;
	}
	private String results,control,pattern,out;
	private long expectedInputs=-1;
	private boolean overwrite;
	private final PrintStream outstream;
	public static final String COUNTS_HEADER="family_id\tassigned_members\tinput_records\teligible_records\ttotal_assigned\ttotal_unassigned\tfiltered_records";
	public static final String REPORT_HEADER="candidate_index\toperation\tfeasible\treason\tinput_a\tinput_b\tinitial_affected_members\tcontrol_affected_members\tmodified_affected_members\tpost_minus_initial\tcontrol_minus_initial\tmodified_minus_matched_control\tinput_records";
}
