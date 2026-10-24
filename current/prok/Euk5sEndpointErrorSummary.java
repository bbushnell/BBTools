package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.IntLongHashMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import structures.ByteBuilder;

/** Read-only reduction of existing one-to-one euk5S grades. Positive signed
 * error is downstream in transcript space; only matched CALL rows enter the
 * distributions. Task26 is a separate historical stratum, never pooled twice.
 * @author Ganyu
 */
public final class Euk5sEndpointErrorSummary {
	public static void main(String[] args) throws Exception {
		String runs=null,out=null,task26=null;int expected=-1;
		for(String arg:args){final int e=arg.indexOf('=');require(e>0,"Expected key=value");
			final String k=arg.substring(0,e).toLowerCase(java.util.Locale.ROOT),v=arg.substring(e+1);
			if(k.equals("runs")){runs=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("task26pairs")){task26=v;}else if(k.equals("expectedgenomes")){expected=Integer.parseInt(v);}
			else if(k.equals("selftest") && v.equals("t")){selfTest();return;}else{throw new IllegalArgumentException(k);}}
		require(runs!=null && out!=null && expected>0,"Required runs, outprefix, expectedgenomes");
		final Euk5sEndpointErrorSummary tool=new Euk5sEndpointErrorSummary();
		tool.loadRuns(runs);require(tool.genomes==expected,"Genome count mismatch: "+tool.genomes+" != "+expected);
		require(tool.groups.size()==7,"Expected six clades and pooled, observed "+tool.groups.size());
		if(task26!=null){tool.loadHistorical(task26);}
		tool.write(out);System.err.println("EUK5S_ENDPOINT_SUMMARY_PASS genomes="+tool.genomes+" groups="+tool.groups.size());
	}

	void loadRuns(String paths){
		final LineParser1 p=new LineParser1('\t');
		for(String path:paths.split(",")){
			final ByteFile in=RrnaSourceTargetPrep.open(path);
			try{header(in,"arm\tclade\tassembly\tgrade\tstages\tcallerSeconds",path);
				for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==6,"Bad RUNS row");
					if(!p.termEquals("production",0) || p.termEquals("prok",1)){continue;}
					final String clade=p.parseString(1),assembly=p.parseString(2),grade=p.parseString(3);
					require(seen.add(clade+"\t"+assembly),"Duplicate clade/assembly: "+clade+"/"+assembly);
					require(grade.endsWith(".summary.tsv"),"Expected grade.summary.tsv");
					final RrnaSeedBenchmarkSummary.Row summary=RrnaSeedBenchmarkSummary.loadGrade(grade,"production",clade,assembly);
					require(summary.family.equals("euk5S") && !summary.negative,"Endpoint cohort must be positive euk5S truth");
					final String loci=grade.substring(0,grade.length()-12)+".loci.tsv.gz";
					final Group g=group(clade),pool=group("pooled");
					loadLoci(loci,assembly,summary,g,pool);g.genomes++;pool.genomes++;genomes++;
					inputs.add(clade+"\t"+assembly+"\t"+grade+"\t"+loci);
				}
			}finally{require(!in.close(),"RUNS read failed: "+path);}
		}
	}

	void loadLoci(String path,String assembly,RrnaSeedBenchmarkSummary.Row summary,Group g,Group pool){
		final ByteFile in=RrnaSourceTargetPrep.open(path);final LineParser1 p=new LineParser1('\t');
		long tp=0,fp=0,truthTp=0,fn=0,e5=0,e3=0,both=0;
		final IntLongHashMap t5=new IntLongHashMap(),t3=new IntLongHashMap(),c5=new IntLongHashMap(),c3=new IntLongHashMap();
		try{header(in,"assembly\tseqid\trole\tindex\tstart\tstop\tstrand\toutcome\tpartnerIndex\td5\td3\tcontrolOverlapAnyStrand",path);
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==12 && p.termEquals(assembly,0),"Locus identity/schema mismatch");
				require(p.termEquals("+",6) || p.termEquals("-",6),"Known strand required");
				final boolean call=p.termEquals("CALL",2),truth=p.termEquals("TRUTH",2),matched=p.termEquals("TP",7);
				require(call || truth,"Unexpected locus role");
				if(matched){require(p.parseInt(8)>=0,"Matched locus lacks partner");final int d5=p.parseInt(9),d3=p.parseInt(10);
					if(call){tp++;c5.increment(d5,1);c3.increment(d3,1);g.add(d5,d3);pool.add(d5,d3);if(d5==0){e5++;}if(d3==0){e3++;}if(d5==0 && d3==0){both++;}}
					else{truthTp++;t5.increment(d5,1);t3.increment(d3,1);}
				}else{require(p.parseInt(8)==-1 && p.termEquals("NA",9) && p.termEquals("NA",10),"Unmatched locus has endpoint/partner data");
					if(call){require(p.termEquals("FP",7),"Unmatched CALL must be FP");fp++;}
					else{require(p.termEquals("FN",7),"Unmatched TRUTH must be FN");fn++;}}
			}
		}finally{require(!in.close(),"Locus read failed: "+path);}
		require(tp==truthTp && tp==summary.tp && fp==summary.fp && fn==summary.fn,"Locus/summary counts disagree: "+path);
		require(e5==summary.exact5 && e3==summary.exact3 && both==summary.both,"Locus/summary exact endpoint counts disagree: "+path);
		require(sameHistogram(c5,t5) && sameHistogram(c3,t3),"CALL/TRUTH paired error histograms disagree: "+path);
	}

	void loadHistorical(String path){
		final ByteFile in=RrnaSourceTargetPrep.open(path);final LineParser1 p=new LineParser1('\t');final Group g=group("task26_fungi112_separate");
		try{header(in,"family\ttruthId\tcallId\tseqid\tstrand\ttruthStart\ttruthStop\tcallStart\tcallStop\texact5\texact3\tbothExact\tdelta5\tdelta3\tabsDelta5\tabsDelta3",path);
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==16,"Bad historical pair schema");if(!p.termEquals("5S",0)){continue;}
				final boolean minus=p.termEquals("-",4);require(minus || p.termEquals("+",4),"Historical strand unknown");
				final int d5=minus?p.parseInt(6)-p.parseInt(8):p.parseInt(7)-p.parseInt(5);
				final int d3=minus?p.parseInt(5)-p.parseInt(7):p.parseInt(8)-p.parseInt(6);
				require(d5==p.parseInt(12) && d3==p.parseInt(13) && Math.abs(d5)==p.parseInt(14) && Math.abs(d3)==p.parseInt(15),"Historical signed errors disagree with coordinates");g.add(d5,d3);
			}
		}finally{require(!in.close(),"Historical pair read failed");}
		g.genomes=112;inputs.add("task26_fungi112_separate\tNA\tNA\t"+path);
	}

	Group group(String name){Group g=byName.get(name);if(g==null){g=new Group(name);byName.put(name,g);groups.add(g);}return g;}
	void write(String prefix){
		RrnaCmWindowParity.fresh(prefix,new String[]{".summary.tsv",".signed.tsv",".inputs.tsv"});
		final ByteStreamWriter summary=RrnaSourceTargetPrep.writer(prefix+".summary.tsv"),hist=RrnaSourceTargetPrep.writer(prefix+".signed.tsv"),input=RrnaSourceTargetPrep.writer(prefix+".inputs.tsv");
		final ByteBuilder b=new ByteBuilder();
		try{summary.println("clade\tgenomes\tend\tmatchedTP\texact\tsignedMean\tabsP95\tabsP98");hist.println("clade\tend\tsignedError\tcount");input.println("clade\tassembly\tgrade\tloci_or_pairs");
			for(String line:inputs){input.println(line);}
			for(Group g:groups){for(int end=0;end<2;end++){
				final IntLongHashMap signed=g.signed[end],absolute=new IntLongHashMap();final int[] keys=signed.toArray();Arrays.sort(keys);long n=0,sum=0;
				for(int d:keys){final long count=signed.get(d);require(d!=Integer.MIN_VALUE && count>0,"Endpoint histogram invalid");n+=count;sum+=count*d;absolute.increment(Math.abs(d),count);
					b.clear().append(g.name).tab().append(end==0?"5prime":"3prime").tab().append(d).tab().append(count).nl();hist.print(b);}
				require(n==g.n,"Histogram must conserve matched CALL count");
				b.clear().append(g.name).tab().append(g.genomes).tab().append(end==0?"5prime":"3prime").tab().append(n).tab().append(signed.get(0)<0?0:signed.get(0)).tab();
				if(n==0){b.append("NA\tNA\tNA");}else{b.append(sum/(double)n,8).tab().append(quantile(absolute,n,95)).tab().append(quantile(absolute,n,98));}b.nl();summary.print(b);
			}}
		}finally{RrnaCmWindowParity.finish(summary,hist,input);}
	}

	static int quantile(IntLongHashMap histogram,long n,int percent){
		require(n>0 && percent>0 && percent<=100,"Nearest-rank quantile requires positive population and percentile");
		final long rank=n/100*percent+(n%100*percent+99)/100;long cumulative=0;final int[] keys=histogram.toArray();Arrays.sort(keys);
		for(int key:keys){cumulative+=histogram.get(key);if(cumulative>=rank){return key;}}
		throw new IllegalStateException("Histogram smaller than declared quantile population");
	}
	static boolean sameHistogram(IntLongHashMap a,IntLongHashMap b){final int[] ak=a.toArray(),bk=b.toArray();if(ak.length!=bk.length){return false;}for(int k:ak){if(a.get(k)!=b.get(k)){return false;}}return true;}
	static void header(ByteFile in,String expected,String path){final byte[] h=in.nextLine();require(h!=null && new String(h,StandardCharsets.UTF_8).equals(expected),"Unexpected header: "+path);}
	static void require(boolean condition,String message){if(!condition){throw new IllegalArgumentException(message);}}
	static void selfTest() throws Exception {
		final IntLongHashMap h=new IntLongHashMap();h.increment(0,96);h.increment(3,4);
		require(quantile(h,100,95)==0 && quantile(h,100,98)==3,"Nearest-rank tail fixture failed");
		final java.nio.file.Path dir=java.nio.file.Files.createTempDirectory("euk5s-endpoint-test-");
		final java.nio.file.Path loci=dir.resolve("loci.tsv"),pairs=dir.resolve("pairs.tsv");
		try{
			final String rows="assembly\tseqid\trole\tindex\tstart\tstop\tstrand\toutcome\tpartnerIndex\td5\td3\tcontrolOverlapAnyStrand\n"
				+"fixture\tplus\tCALL\t0\t102\t120\t+\tTP\t0\t2\t-3\t0\n"
				+"fixture\tplus\tTRUTH\t0\t100\t123\t+\tTP\t0\t2\t-3\t0\n"
				+"fixture\tminus\tCALL\t0\t205\t220\t-\tTP\t0\t2\t-5\t0\n"
				+"fixture\tminus\tTRUTH\t0\t200\t222\t-\tTP\t0\t2\t-5\t0\n"
				+"fixture\tplus\tCALL\t1\t300\t420\t+\tFP\t-1\tNA\tNA\t0\n"
				+"fixture\tminus\tTRUTH\t1\t500\t620\t-\tFN\t-1\tNA\tNA\t0\n";
			java.nio.file.Files.write(loci,rows.getBytes(StandardCharsets.UTF_8));
			final Euk5sEndpointErrorSummary tool=new Euk5sEndpointErrorSummary();
			final RrnaSeedBenchmarkSummary.Row expected=new RrnaSeedBenchmarkSummary.Row();expected.tp=2;expected.fp=1;expected.fn=1;
			final Group g=new Group("fixture"),pool=new Group("pooled");tool.loadLoci(loci.toString(),"fixture",expected,g,pool);
			require(g.n==2 && pool.n==2 && g.signed[0].get(2)==2 && g.signed[1].get(-3)==1 && g.signed[1].get(-5)==1,"Only matched CALL rows enter strand-aware histograms");
			expected.tp=3;boolean rejected=false;
			try{tool.loadLoci(loci.toString(),"fixture",expected,new Group("bad"),new Group("badpool"));}catch(IllegalArgumentException e){rejected=true;}
			require(rejected,"Locus/summary mismatch must fail");
			final String historical="family\ttruthId\tcallId\tseqid\tstrand\ttruthStart\ttruthStop\tcallStart\tcallStop\texact5\texact3\tbothExact\tdelta5\tdelta3\tabsDelta5\tabsDelta3\n"
				+"5S\tT0\tC0\tplus\t+\t100\t123\t102\t120\tfalse\tfalse\tfalse\t2\t-3\t2\t3\n"
				+"5S\tT1\tC1\tminus\t-\t200\t222\t205\t220\tfalse\tfalse\tfalse\t2\t-5\t2\t5\n";
			java.nio.file.Files.write(pairs,historical.getBytes(StandardCharsets.UTF_8));tool.loadHistorical(pairs.toString());
			require(tool.group("task26_fungi112_separate").n==2,"Historical plus/minus pairs must be read");
			java.nio.file.Files.write(pairs,historical.replace("\t2\t-5\t2\t5\n","\t-2\t-5\t2\t5\n").getBytes(StandardCharsets.UTF_8));rejected=false;
			try{new Euk5sEndpointErrorSummary().loadHistorical(pairs.toString());}catch(IllegalArgumentException e){rejected=true;}
			require(rejected,"Wrong minus-strand signed delta must fail");
		}finally{java.nio.file.Files.deleteIfExists(loci);java.nio.file.Files.deleteIfExists(pairs);java.nio.file.Files.deleteIfExists(dir);}
		System.err.println("EUK5S_ENDPOINT_SELFTEST_PASS paired_strands unmatched_exclusion count_mismatch signed_delta_mismatch quantiles");
	}
	public static final class Group{Group(String n_){name=n_;}void add(int d5,int d3){signed[0].increment(d5,1);signed[1].increment(d3,1);n++;}final String name;public long n,genomes;public final IntLongHashMap[] signed={new IntLongHashMap(),new IntLongHashMap()};}
	final ObjectMap<String,Group> byName=new ObjectMap<String,Group>(String.class,Group.class);
	final ObjectSet<String> seen=new ObjectSet<String>(String.class);
	final ArrayList<Group> groups=new ArrayList<Group>();final ArrayList<String> inputs=new ArrayList<String>();long genomes=0;
}
