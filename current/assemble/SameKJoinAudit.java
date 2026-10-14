package assemble;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Audits saved same-K candidate provenance without loading reads or reference.
 * Original walk stop codes and freshly queried endpoint junction flags remain
 * separate: a backward/hidden branch can stop before the next base is appended.
 * @author Fischl
 */
public final class SameKJoinAudit {

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs a single streaming pair audit through native argument and I/O helpers. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SameKJoinAudit.class, false);
		try{new SameKJoinAudit(pp.args, pp.outstream).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires original candidate and label files and one new audit output. */
	private SameKJoinAudit(final String[] args, final PrintStream outstream_){
		outstream=outstream_;
		final Parser parser=new Parser();
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			String name=split[0].toLowerCase(Locale.ROOT);
			while(name.startsWith("-")){name=name.substring(1);}
			final String value=split.length>1 ? split[1] : null;
			if(name.equals("labels")){labels=value;}
			else if(!parser.parse(arg, name, value)){throw new IllegalArgumentException("Unknown audit argument: "+arg);}
		}
		input=parser.in1; output=parser.out1;
		if(input==null || labels==null || output==null || parser.in2!=null || parser.out2!=null){
			throw new IllegalArgumentException("Required in=candidates.tsv labels=labels.tsv out=audit.tsv; no paired FASTQ inputs.");
		}
		TadpoleGraph.checkPaths(input, new ArrayList<String>(Arrays.asList(input, labels)), output);
		if(!Tools.testInputFiles(false, true, input, labels) || !Tools.testOutputFiles(false, false, false, output)){
			throw new IllegalArgumentException("Audit inputs must exist and output must be new.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts source/destination provenance in lockstep, without changing any labels. */
	private void run(){
		assert(input!=null && labels!=null) : "Every census row must retain its original paired truth record.";
		final long[] counts=new long[ShaveObject.MAX_CODE*ShaveObject.MAX_CODE*4*3];
		final ByteFile cf=ByteFile.makeByteFile(input, true), lf=ByteFile.makeByteFile(labels, true);
		final LineParser1 cp=new LineParser1('\t'), lp=new LineParser1('\t');
		int event=0, k=0;
		long recordedBranch=0, neitherFlag=0, branchWithoutFlag=0, ownership=0;
		boolean error=false;
		try{
			requireHeader(cf.nextLine(), SameKJoinDiagnostic.HEADER);
			final byte[] comment=lf.nextLine();
			if(comment==null || comment.length==0 || comment[0]!='#'){throw new IllegalArgumentException("Missing truth schema comment.");}
			requireHeader(lf.nextLine(), LABEL_HEADER);
			for(byte[] row=cf.nextLine(); row!=null; row=cf.nextLine()){
				final byte[] label=lf.nextLine();
				if(label==null){throw new IllegalArgumentException("Missing paired label.");}
				cp.set(row); lp.set(label); event++;
				if(cp.terms()!=36 || lp.terms()!=17 || cp.parseInt(0)!=event || lp.parseInt(0)!=event){
					throw new IllegalArgumentException("Schema/ID mismatch at event "+event);
				}
				final int current=cp.parseInt(1);
				if(current<2 || current!=lp.parseInt(1) || (k!=0 && current!=k)){throw new IllegalArgumentException("Mixed or mismatched K values.");}
				k=current;
				final int sc=cp.parseInt(10), dc=cp.parseInt(11), sj=cp.parseInt(14), dj=cp.parseInt(15);
				if(sc<0 || sc>=ShaveObject.MAX_CODE || dc<0 || dc>=ShaveObject.MAX_CODE || (sj&~1)!=0 || (dj&~1)!=0){
					throw new IllegalArgumentException("Invalid stop code/junction flag at event "+event);
				}
				final int truth=lp.termEquals("supported", 2) ? 0 : lp.termEquals("contradicted", 2) ? 1 : lp.termEquals("unresolved", 2) ? 2 : -1;
				if(truth<0 || !lp.termEquals(TARGETS[truth], 3)){throw new IllegalArgumentException("Unknown or inconsistent truth label/target.");}
				counts[index(sc, dc, sj, dj, truth)]++;
				final boolean branch=ShaveObject.isBranchCode(sc) || ShaveObject.isBranchCode(dc);
				if(branch){recordedBranch++;}
				if(sj==0 && dj==0){neitherFlag++; if(branch){branchWithoutFlag++;}}
				if(sc==ShaveObject.BAD_OWNER || dc==ShaveObject.BAD_OWNER){ownership++;}
			}
			if(lf.nextLine()!=null){throw new IllegalArgumentException("Extra labels after final candidate.");}
		}finally{error=cf.close(); error|=lf.close();}
		if(error){throw new RuntimeException("Error reading audit inputs.");}
		write(counts, k, event);
		outstream.println("SAME_K_AUDIT_PASS candidates="+event+" recorded_branch="+recordedBranch+
				" neither_endpoint_junction="+neitherFlag+" branch_without_endpoint_junction="+branchWithoutFlag+" ownership="+ownership);
	}

	/** Writes sparse observed groups while verifying exact count conservation. */
	private void write(final long[] counts, final int k, final int expected){
		assert(counts.length==ShaveObject.MAX_CODE*ShaveObject.MAX_CODE*12) : "One bin per code pair, flag pair and truth class.";
		final ByteStreamWriter writer=new ByteStreamWriter(output, false, false, false);
		writer.start();
		final ByteBuilder row=new ByteBuilder();
		long total=0;
		try{
			writer.print("k\tsource_code\tdest_code\tsource_junction\tdest_junction\tlabel\tcount\n");
			for(int sc=0; sc<ShaveObject.MAX_CODE; sc++){
				for(int dc=0; dc<ShaveObject.MAX_CODE; dc++){
					for(int sj=0; sj<2; sj++){
						for(int dj=0; dj<2; dj++){
							for(int label=0; label<3; label++){
								final long n=counts[index(sc, dc, sj, dj, label)];
								if(n==0){continue;}
								total+=n;
								writer.print(row.clear().append(k).tab().append(sc).tab().append(dc).tab().append(sj).tab()
										.append(dj).tab().append(LABELS[label]).tab().append(n).nl());
							}
						}
					}
				}
			}
			if(total!=expected){throw new IllegalStateException("Audit grouping lost records: "+total+" != "+expected);}
		}finally{if(writer.poisonAndWait()){throw new RuntimeException("Audit output failed.");}}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Fixed-width primitive bins avoid maps or one object per candidate. */
	private static int index(final int sc, final int dc, final int sj, final int dj, final int label){
		assert(sc>=0 && dc>=0 && sc<ShaveObject.MAX_CODE && dc<ShaveObject.MAX_CODE
				&& sj>=0 && sj<2 && dj>=0 && dj<2 && label>=0 && label<3) : "Validated provenance must fit its fixed schema.";
		return ((((sc*ShaveObject.MAX_CODE+dc)*2+sj)*2+dj)*3+label);
	}

	/** A shifted truth or candidate column is an error, never a new control class. */
	private static void requireHeader(final byte[] row, final String expected){
		if(row==null || !expected.equals(new String(row, StandardCharsets.US_ASCII))){throw new IllegalArgumentException("Unexpected audit input schema.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final PrintStream outstream;
	private String input, labels, output;
	private static final String[] LABELS={"supported", "contradicted", "unresolved"}, TARGETS={"1", "0", "NA"};
	private static final String LABEL_HEADER="event\tk\tlabel\ttarget\twhole_hits\tleft_hits\tright_hits\tsource_hits\tsource_unknown\tsource_left_mask\tsource_right_mask\tsource_pair_mask\tdest_hits\tdest_unknown\tdest_left_mask\tdest_right_mask\tdest_pair_mask";
}
