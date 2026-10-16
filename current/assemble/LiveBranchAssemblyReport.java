package assemble;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Summarizes two named QUAST arms, retaining empty and undefined outcomes. @author Fischl */
public final class LiveBranchAssemblyReport {

	/** Uses standard config expansion and rejects unknown flags. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, LiveBranchAssemblyReport.class, false);
		try{new LiveBranchAssemblyReport(pp.args, pp.outstream).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires a new output directory and an explicit expected panel size. */
	private LiveBranchAssemblyReport(final String[] args, final PrintStream outstream_){
		outstream=outstream_;
		final Parser parser=new Parser();
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			String key=split[0].toLowerCase(Locale.ROOT);
			while(key.startsWith("-")){key=key.substring(1);}
			final String value=split.length>1 ? split[1] : null;
			if(key.equals("expected")){expected=Integer.parseInt(value);}
			else if(key.equals("arms")){arms=value==null ? new String[0] : value.split(",", -1);}
			else if(!parser.parse(arg, key, value)){throw new IllegalArgumentException("Unknown report argument: "+arg);}
		}
		input=parser.in1; output=parser.out1;
		if(arms.length!=2 || !arms[0].matches("[a-z][a-z0-9_]*") ||
			!arms[1].matches("[a-z][a-z0-9_]*") || arms[0].equals(arms[1])){
			throw new IllegalArgumentException("arms requires two distinct lowercase labels; first is the baseline.");
		}
		if(input==null || output==null || expected<1 || !new File(input).isDirectory() || new File(output).exists()){
			throw new IllegalArgumentException("Required in=collected_directory out=NEW_directory expected=genomes.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Validates the complete frozen VALID panel before creating output. */
	private void run(){
		final ArrayList<String> ids=new ArrayList<String>();
		final ArrayList<double[][]> panel=new ArrayList<double[][]>();
		final HashSet<String> seen=new HashSet<String>();
		final ByteFile manifest=ByteFile.makeByteFile(input+"/manifest.tsv", true);
		final LineParser1 lp=new LineParser1('\t');
		try{
			lp.set(manifest.nextLine());
			if(lp.terms()!=7 || !lp.termEquals("genome", 0) || !lp.termEquals("fold", 1) || !lp.termEquals("run", 6)){
				throw new IllegalArgumentException("Expected the seven-column live-validation manifest.");
			}
			for(byte[] line=manifest.nextLine(); line!=null; line=manifest.nextLine()){
				lp.set(line);
				if(lp.terms()!=7 || !lp.termEquals("validation", 1) || lp.parseInt(2)<24 || lp.parseInt(2)>127 || lp.parseInt(3)<1){
					throw new IllegalArgumentException("Wrong fold, K, depth or manifest width.");
				}
				final String id=lp.parseString(0);
				if(!id.matches("g[0-9]{4}") || !seen.add(id)){throw new IllegalArgumentException("Unsafe/duplicate genome: "+id);}
				ids.add(id); panel.add(readReport(input+"/runs/"+id, arms));
			}
		}finally{if(manifest.close()){throw new IllegalStateException("Manifest read failed.");}}
		if(ids.size()!=expected){throw new IllegalArgumentException("Expected "+expected+" genomes, found "+ids.size());}
		if(!new File(output).mkdir()){throw new IllegalArgumentException("Cannot create output: "+output);}
		writeRows(ids, panel); writeSummary(panel);
		outstream.println("LIVE_BRANCH_ASSEMBLY_REPORT_PASS genomes="+ids.size()+" arms=2");
	}

	/** Emits each genome/arm, including explicit undefined-NGA50 flags. */
	private void writeRows(final ArrayList<String> ids, final ArrayList<double[][]> panel){
		assert(ids.size()==panel.size()) : "Each validated manifest genome supplies exactly one paired report.";
		final ByteStreamWriter out=open(output+"/metrics.tsv");
		final ByteBuilder row=new ByteBuilder();
		try{
			out.print("genome\tarm\tcontigs\tnga50_zero\tglobal\tlocal\tgenome_fraction\tbases\tnga50_defined\n");
			for(int g=0; g<ids.size(); g++){
				for(int a=0; a<2; a++){
					row.clear().append(ids.get(g)).tab().append(arms[a]);
					for(double value : panel.get(g)[a]){row.tab().append(value, 6);}
					out.print(row.nl());
				}
			}
		}finally{if(out.poisonAndWait()){throw new IllegalStateException("Metrics write failed.");}}
	}

	/** Every genome has equal weight; missing NGA50 contributes zero, never a dropped row. */
	private void writeSummary(final ArrayList<double[][]> panel){
		assert(!panel.isEmpty()) : "A complete nonempty manifest is required for means and medians.";
		final ByteStreamWriter out=open(output+"/summary.tsv");
		final ByteBuilder row=new ByteBuilder();
		try{
			out.print("arm\tgenomes\tundefined_nga50\tglobal_total\tlocal_total\tcontigs_total\tcontigs_median\tnga50_zero_mean\tnga50_zero_median\tgenome_fraction_mean\tnga50_better_"+arms[0]+
				"\tnga50_equal_"+arms[0]+"\tnga50_worse_"+arms[0]+"\tglobal_more_"+arms[0]+"\tglobal_less_"+arms[0]+"\n");
			for(int a=0; a<2; a++){
				final double[] nga=new double[panel.size()], counts=new double[panel.size()];
				long global=0, local=0, contigs=0;
				int undefined=0, better=0, equal=0, worse=0, moreErrors=0, fewerErrors=0;
				double sum=0, gf=0;
				for(int g=0; g<panel.size(); g++){
					final double[] values=panel.get(g)[a], baseline=panel.get(g)[0];
					counts[g]=values[0]; nga[g]=values[1];
					contigs+=(long)values[0]; global+=(long)values[2]; local+=(long)values[3];
					sum+=values[1]; gf+=values[4]; undefined+=values[6]==0 ? 1 : 0;
					if(values[1]>baseline[1]){better++;}else if(values[1]<baseline[1]){worse++;}else{equal++;}
					if(values[2]>baseline[2]){moreErrors++;}else if(values[2]<baseline[2]){fewerErrors++;}
				}
				row.clear().append(arms[a]).tab().append(panel.size()).tab().append(undefined).tab().append(global).tab()
					.append(local).tab().append(contigs).tab().append(SameKAssemblyReport.median(counts), 6).tab()
					.append(sum/panel.size(), 6).tab().append(SameKAssemblyReport.median(nga), 6).tab().append(gf/panel.size(), 6)
					.tab().append(better).tab().append(equal).tab().append(worse).tab().append(moreErrors).tab().append(fewerErrors);
				out.print(row.nl());
			}
		}finally{if(out.poisonAndWait()){throw new IllegalStateException("Summary write failed.");}}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Empty assemblies require explicit markers and must be absent from QUAST's header. */
	private static double[][] readReport(final String directory, final String[] arms){
		assert(directory!=null) : "The report directory is constructed from a validated genome ID.";
		assert(arms.length==2) : "Two validated arm labels map QUAST columns and explicit empty markers.";
		final boolean[] empty={new File(directory+"/"+arms[0]+".empty").isFile(), new File(directory+"/"+arms[1]+".empty").isFile()};
		final double[][] data=new double[2][7];
		final int n=(empty[0] ? 0 : 1)+(empty[1] ? 0 : 1);
		final String path=directory+"/quast/report.tsv";
		if(n==0){
			if(new File(path).exists()){throw new IllegalArgumentException("QUAST report contradicts two empty markers: "+directory);}
			return data;
		}
		final ByteFile file=ByteFile.makeByteFile(path, true);
		final LineParser1 lp=new LineParser1('\t');
		final boolean[] found=new boolean[METRICS.length];
		try{
			lp.set(file.nextLine());
			if(lp.terms()!=n+1 || !lp.termEquals("Assembly", 0)){throw new IllegalArgumentException("Wrong QUAST header: "+path);}
			for(int a=0, column=1; a<2; a++){
				if(!empty[a] && !lp.termEquals(arms[a], column++)){throw new IllegalArgumentException("Wrong QUAST arm order: "+path);}
			}
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				lp.set(line);
				if(lp.terms()!=n+1){throw new IllegalArgumentException("Wrong QUAST row width: "+path);}
				for(int m=0; m<METRICS.length; m++){
					if(!lp.termEquals(METRICS[m], 0)){continue;}
					if(found[m]){throw new IllegalArgumentException("Duplicate metric: "+METRICS[m]);}
					found[m]=true;
					for(int a=0, column=1; a<2; a++){
						if(empty[a]){continue;}
						final boolean missing=lp.termEquals("-", column);
						if(missing && m!=1){throw new IllegalArgumentException("Undefined required metric: "+METRICS[m]);}
						final double value=missing ? 0 : lp.parseDouble(column);
						if(!Double.isFinite(value) || value<0 || (m!=4 && value!=Math.rint(value)) || (m==4 && value>100)){
							throw new IllegalArgumentException("Invalid QUAST metric: "+METRICS[m]+"="+value);
						}
						data[a][m]=value;
						if(m==1){data[a][6]=missing ? 0 : 1;}
						column++;
					}
				}
			}
		}finally{if(file.close()){throw new IllegalStateException("QUAST read failed.");}}
		for(int m=0; m<found.length; m++){if(!found[m]){throw new IllegalArgumentException("Missing metric: "+METRICS[m]);}}
		return data;
	}

	/** Output files are new and every caller checks the writer's terminal error state. */
	private static ByteStreamWriter open(final String path){
		final ByteStreamWriter out=new ByteStreamWriter(path, false, false, false);
		out.start(); return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final String input, output;
	private final PrintStream outstream;
	private int expected=-1;
	private String[] arms={"native", "live"};
	private static final String[] METRICS={"# contigs", "NGA50", "# misassemblies", "# local misassemblies", "Genome fraction (%)", "Total length"};
}
