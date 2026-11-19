package prot.cluster;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntHashMap;

/** Descriptive in-family tail-GC retention from the accepted SHORT annotation.
 * No DNA inference, panel scoring, clustering or provenance hashing.
 * @author Brian, Qiqi
 */
public final class FamilyTailGcSweep{

	public static void main(String[] args){
		try{run(args);}catch(Throwable failure){failure.printStackTrace();System.exit(1);}
	}
	public static void run(String[] args)throws Exception{
		final PreParser pp=new PreParser(args,FamilyTailGcSweep.class,false);
		final int priorThreads=Shared.threads(),priorBuffers=Shared.numBuffers();
		// The native reader selects inline gzip inflation at one shared thread.
		// Bind that budget even when called from a larger fixture JVM, then restore it.
		Shared.setThreads(1);
		try{new FamilyTailGcSweep(pp.args).process();}
		finally{Shared.setThreads(priorThreads);Shared.setBuffers(priorBuffers);Shared.closeStream(pp.outstream);}
	}
	private FamilyTailGcSweep(String[] args){
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			require(pair.length==2 && !pair[1].isEmpty(),"Expected flag=value");
			require(options.put(pair[0].toLowerCase(Locale.ROOT),pair[1])==null,"Duplicate argument");
		}
		queryPath=take(options,"querygc");familyPath=take(options,"families");memberPath=take(options,"members");
		out=Paths.get(take(options,"out")).toAbsolutePath().normalize();
		expectedQueries=positive(take(options,"expectedqueries"));
		expectedAssigned=positive(take(options,"expectedassigned"));
		expectedFamilies=Math.toIntExact(positive(take(options,"expectedfamilies")));
		require(take(options,"multiples").equals(GRID),"The reviewed finite sweep is required");
		require(take(options,"verification").equals("structural"),"Explicit inherited/structural provenance required");
		require(options.isEmpty(),"Unknown arguments: "+options.keySet());
		require(expectedAssigned<=expectedQueries && expectedFamilies<=expectedAssigned,"Invalid expected cohort");
		require(!Files.exists(out,LinkOption.NOFOLLOW_LINKS),"Output must be a fresh directory");
	}

	private void process()throws Exception{
		for(String path:new String[]{queryPath,familyPath,memberPath}){inputs.add(new Snapshot(Paths.get(path)));}
		readFamilies();Files.createDirectory(out);
		final long wall1=System.nanoTime(),cpu1=cpuTime();
		scan(true);final long cpu2=cpuTime(),wall2=System.nanoTime();
		reconcile();checkInputs();
		final long wall3=System.nanoTime(),cpu3=cpuTime();
		scan(false);final long cpu4=cpuTime(),wall4=System.nanoTime();
		for(Family f:families){
			require(f.visited==f.n,"Second-pass family population changed");
			for(int i=0;i<MULTIPLES.length;i++){
				require(f.excluded[i]>=0 && f.excluded[i]<=f.n,"Invalid exclusion count");
				require(i==0 || f.excluded[i]<=f.excluded[i-1],"Gate acceptance must be nested");
			}
		}
		checkInputs();writeReports();
		writeText("phase_times.tsv","phase\twall_seconds\tjava_cpu_seconds\n"+
			"moments\t"+seconds(wall2-wall1)+"\t"+seconds(cpu2-cpu1)+"\n"+
			"decisions\t"+seconds(wall4-wall3)+"\t"+seconds(cpu4-cpu3)+"\n");
		final ByteStreamWriter snapshots=writer("inputs.tsv");
		try{
			snapshots.println("path\tbytes\tmtime\tfile_key\tverification");
			for(Snapshot s:inputs){snapshots.println(s.path+"\t"+s.before.size()+"\t"+s.before.lastModifiedTime()+"\t"+
				String.valueOf(s.before.fileKey())+"\tmetadata_stability_inherited_content_identity");}
		}finally{require(!snapshots.poisonAndWait(),"Input receipt write failed");}
		checkInputs();
		writeText("evidence.paths",out.resolve("family_gate.tsv")+"\n"+out.resolve("summary.tsv")+"\n"+
			out.resolve("coverage.tsv")+"\n"+out.resolve("phase_times.tsv")+"\n"+out.resolve("inputs.tsv")+"\n"+
			out.resolve("evidence.paths")+"\n"+out.resolve("PASS")+"\n");
		writeText("PASS","FAMILY_TAIL_GC_SWEEP_PASS\n");
		System.out.println("FAMILY_TAIL_GC_SWEEP_PASS queries="+expectedQueries+" assigned="+expectedAssigned+" families="+families.size());
	}

	/** Recover only a unique integer in a conservative one-unit decimal interval.
	 * A separate exact producer-formatter roundtrip is mandatory for each row.
	 * All arithmetic that could overflow fails; no rounded-fraction fallback.
	 */
	public static long recoverNumerator(long decimalScaled,long denominator){
		require(decimalScaled>=0 && decimalScaled<=SCALE && denominator>0,"Invalid measured GC fraction/count");
		require(Math.multiplyExact(denominator,2)<SCALE,"Nonunique numerator interval");
		final long product=Math.multiplyExact(decimalScaled,denominator);
		final long candidate=Math.addExact(product,SCALE/2)/SCALE;
		require(candidate>=0 && candidate<=denominator,"Numerator outside canonical count");
		final long distance=Math.abs(Math.subtractExact(product,Math.multiplyExact(candidate,SCALE)));
		require(distance<=denominator,"No integer numerator in decimal interval");
		return candidate;
	}

	private void readFamilies(){
		final ByteFile bf=ByteFile.makeByteFile1(familyPath,false);
		final LineParser1 p=new LineParser1('\t');final HashSet<String> labels=new HashSet<String>();
		long assigned=0;
		try{
			header(bf,FAMILY_HEADER);
			for(byte[] row=bf.nextLine();row!=null;row=bf.nextLine()){
				p.set(row);require(p.terms()==21,"Family report requires21 fields");
				final int id=p.parseInt(0),rank=p.parseInt(1);final String label=p.parseString(2);
				final long n=nonnegative(p,3);
				require(id>=0 && rank==families.size() && n>0 && byId.get(id)<0 && !label.isEmpty() && labels.add(label),
					"Duplicate/empty family or nonconsecutive roster");
				require(nonnegative(p,4)==n && nonnegative(p,5)==n && nonnegative(p,6)==n &&
					nonnegative(p,7)==0 && nonnegative(p,8)==0 && nonnegative(p,9)==0,"Family GC coverage must be complete");
				final long canonical=Math.addExact(Math.addExact(nonnegative(p,11),nonnegative(p,12)),
					Math.addExact(nonnegative(p,13),nonnegative(p,14)));
				final long gc=Math.addExact(nonnegative(p,12),nonnegative(p,13));
				final Family f=new Family(id,rank,label,n,canonical,gc,nonnegative(p,15),nonnegative(p,10));
				f.expectedMean=scaled(p,16);f.expectedMin=scaled(p,19);f.expectedMax=scaled(p,20);
				if(n==1){require(p.termEquals("NA",17) && p.termEquals("NA",18),"Singleton variance must be unavailable");}
				else{f.expectedVariance=scaled(p,17);f.expectedSd=scaled(p,18);}
				for(int field=16;field<=20;field++){if(p.length(field)==1){f.integerForms|=1<<(field-16);}}
				require(canonical>0 && f.expectedAmbiguousGenes<=n,"Invalid family base/ambiguity totals");
				byId.put(id,rank);families.add(f);assigned=Math.addExact(assigned,n);
			}
		}finally{require(!bf.close(),"Family report read failed");}
		require(families.size()==expectedFamilies && assigned==expectedAssigned,"Family roster cohort differs");
	}

	/** The first pass exact-joins original membership; the second reuses guarded annotations. */
	private void scan(boolean moments){
		final ByteFile qf=ByteFile.makeByteFile1(queryPath,false);
		final ByteFile mf=moments ? ByteFile.makeByteFile1(memberPath,false) : null;
		final LineParser1 p=new LineParser1('\t'),m=new LineParser1('\t');
		final ByteBuilder formatted=new ByteBuilder();
		long rows=0,assigned=0,unassigned=0,unknownUnassigned=0;
		try{
			header(qf,QUERY_HEADER);if(mf!=null){header(mf,MEMBER_HEADER);}
			for(byte[] row=qf.nextLine();row!=null;row=qf.nextLine()){
				p.set(row);require(p.terms()==14 && p.parseLong(0)==rows && p.length(1)>0,"Invalid annotation index/schema");
				if(mf!=null){
					final byte[] member=mf.nextLine();require(member!=null,"Missing original membership row");m.set(member);
					require(m.terms()==4 && m.parseLong(0)==rows && equal(p,1,m,1) && equal(p,2,m,2) && equal(p,3,m,3),
						"Annotation differs from exact original membership");
				}
				final int id=p.parseInt(2);final boolean isAssigned=p.termEquals("assigned",3);
				require(isAssigned || p.termEquals("unassigned",3),"Unexpected filtered/unknown membership status");
				require(isAssigned ? id>=0 : id==-1,"Family ID/status mismatch");
				require(equal(p,4,p,5),"SHORT gene and tail fractions must be identical");
				final boolean measured=!p.termEquals("NA",4);
				long canonical=0,ambiguous=0,gc=0;
				if(measured){
					canonical=nonnegative(p,8);ambiguous=nonnegative(p,9);
					gc=recoverNumerator(scaled(p,4),canonical);
					formatted.clear().append(gc/(double)canonical,12);
					require(matches(p,4,formatted),"Recovered numerator does not roundtrip through original formatter");
				}else{
					require(!isAssigned,"Assigned GC is unexpectedly unavailable");
					if(p.termEquals("NA",8)){require(p.termEquals("NA",9),"Missing canonical/ambiguity counts must agree");}
					else{require(nonnegative(p,8)==0,"Unavailable GC with positive canonical count");nonnegative(p,9);}
					unknownUnassigned++;
				}
				if(isAssigned){
					final int rank=byId.get(id);require(rank>=0,"Assigned family absent from roster");
					final Family f=families.get(rank);final double x=gc/(double)canonical;
					if(moments){
						f.n=Math.addExact(f.n,1);final double delta=x-f.mean;f.mean+=delta/f.n;f.m2+=delta*(x-f.mean);
						f.min=Math.min(f.min,x);f.max=Math.max(f.max,x);
						f.canonical=Math.addExact(f.canonical,canonical);f.gc=Math.addExact(f.gc,gc);
						f.ambiguous=Math.addExact(f.ambiguous,ambiguous);if(ambiguous>0){f.ambiguousGenes++;}
					}else{
						f.visited++;
						if(f.n>=2){for(int i=0;i<MULTIPLES.length;i++){if(Math.abs(x-f.mean)>MULTIPLES[i]*f.sd){f.excluded[i]++;}}}
					}
					assigned++;
				}else{unassigned++;}
				rows++;
				if(rows%1000000==0){System.err.println("TAIL_GC_SWEEP_PROGRESS pass="+(moments ? 1 : 2)+" queries="+rows+" assigned="+assigned);}
			}
			require(mf==null || mf.nextLine()==null,"Extra original membership rows");
		}finally{
			final boolean qe=qf.close(),me=mf!=null && mf.close();require(!qe && !me,"Input stream close failed");
		}
		require(rows==expectedQueries && assigned==expectedAssigned && unassigned==expectedQueries-expectedAssigned,"Query cohort differs");
		if(moments){unassignedUnknown=unknownUnassigned;}else{require(unassignedUnknown==unknownUnassigned,"Neutral coverage changed");}
	}

	private void reconcile(){
		final ByteBuilder text=new ByteBuilder();
		for(Family f:families){
			require(f.n==f.expected && f.canonical==f.expectedCanonical && f.gc==f.expectedGc &&
				f.ambiguous==f.expectedAmbiguous && f.ambiguousGenes==f.expectedAmbiguousGenes,"Family counts/base totals do not reconcile");
			require(Double.isFinite(f.m2) && f.m2>=0,"Invalid centered family moment");
			require(sameRendered(f.mean,f.expectedMean,f.integerForms,0,text) &&
				sameRendered(f.min,f.expectedMin,f.integerForms,3,text) &&
				sameRendered(f.max,f.expectedMax,f.integerForms,4,text),"Family mean/range does not reproduce retained report");
			if(f.n>=2){
				final double variance=f.m2/(f.n-1);f.sd=Math.sqrt(variance);
				require(sameRendered(variance,f.expectedVariance,f.integerForms,1,text) &&
					sameRendered(f.sd,f.expectedSd,f.integerForms,2,text),"Family sample variance/SD differs");
			}
		}
	}
	private void writeReports(){
		final ByteStreamWriter fw=writer("family_gate.tsv"),sw=writer("summary.tsv");
		final ByteBuilder b=new ByteBuilder();
		// all plus six mutually exclusive predeclared size strata; last grid slot is no gate.
		final long[][][] totals=new long[7][MULTIPLES.length+1][8];
		try{
			b.append("stable_family_id\troster_rank\tfamily_id\tmembers\tmean_available\tmean_gc\tsample_sd\tgate_available\tgate_status");
			for(String n:GRID.split(",")){b.append("\texcluded_N").append(n);}fw.println(b);
			for(Family f:families){
				b.clear().append(f.id).tab().append(f.rank).tab().append(f.label).tab().append(f.n)
					.append("\ttrue\t").append(Double.toString(f.mean)).tab();
				if(f.n==1){b.append("NA");}else{b.append(Double.toString(f.sd));}
				b.tab().append(f.n>=2).tab().append(f.n==1 ? "singleton_neutral" : f.sd==0 ? "zero_sd" : "available");
				for(long count:f.excluded){b.tab().append(count);}fw.println(b);
				for(int group:new int[]{0,sizeGroup(f.n)}){
					for(int i=0;i<=MULTIPLES.length;i++){
						final long excluded=i==MULTIPLES.length ? 0 : f.excluded[i];final long[] t=totals[group][i];
						t[0]+=f.n;if(f.n>=2){t[1]+=f.n;}else{t[2]+=f.n;}
						if(f.n>=2 && f.sd==0){t[3]+=f.n;}t[4]+=excluded;t[5]++;if(excluded>0){t[6]++;}
						if(f.n>=2 && f.sd==0){t[7]++;}
					}
				}
			}
			sw.println("size_stratum\tgate\tN\tassigned\tgate_eligible\tsingleton_neutral\tzero_sd_members\texcluded\tretained\texcluded_fraction_all\texcluded_fraction_eligible\tfamilies\tfamilies_with_exclusions\tzero_sd_families");
			for(int group=0;group<totals.length;group++){
				for(int i=0;i<=MULTIPLES.length;i++){
					final long[] t=totals[group][i];require(t[0]==t[1]+t[2] && t[4]<=t[1],"Summary count conservation failed");
					b.clear().append(SIZES[group]).tab().append(i==MULTIPLES.length ? "off\tNA" : "on\t"+MULTIPLES[i]);
					for(int j=0;j<5;j++){b.tab().append(t[j]);}b.tab().append(t[0]-t[4]);
					ratio(b,t[4],t[0]);ratio(b,t[4],t[1]);b.tab().append(t[5]).tab().append(t[6]).tab().append(t[7]);sw.println(b);
				}
			}
			require(totals[0][0][0]==expectedAssigned,"Summary population differs");
		}finally{
			final boolean fe=fw.poisonAndWait(),se=sw.poisonAndWait();require(!fe && !se,"Report output failed");
		}
		writeText("coverage.tsv","metric\tvalue\nqueries\t"+expectedQueries+"\nassigned\t"+expectedAssigned+
			"\nunassigned\t"+(expectedQueries-expectedAssigned)+"\nunassigned_gc_unavailable\t"+unassignedUnknown+
			"\nassigned_gc_unavailable\t0\nfamilies\t"+families.size()+
			"\ninterpretation\tin_family_descriptive_retention_not_held_out_detection\n"+
			"panel_scoring\tnot_performed\nboundary\tclosed_no_epsilon\nsingleton\tneutral_sample_sd_unavailable\n"+
			"zero_sd\texact_mean_only\nfraction_source\tunique_integer_recovery_and_original_formatter_roundtrip\n");
	}

	private static int sizeGroup(long n){return n==1 ? 1 : n<5 ? 2 : n<10 ? 3 : n<50 ? 4 : n<100 ? 5 : 6;}
	private static void ratio(ByteBuilder b,long numerator,long denominator){b.tab();if(denominator==0){b.append("NA");}else{b.append(numerator/(double)denominator,12);}}
	private static long scaled(LineParser1 p,int field){
		if(p.length(field)==1){
			final int digit=p.parseByte(field,0)-'0';require(digit==0 || digit==1,"Canonical integral GC/moment must be0 or1");
			return digit*SCALE;
		}
		require(p.length(field)==14 && (p.parseByte(field,0)=='0' || p.parseByte(field,0)=='1') &&
			p.parseByte(field,1)=='.',"Expected canonical integer or twelve-decimal GC/moment");
		long value=p.parseByte(field,0)-'0';
		for(int i=2;i<14;i++){final int digit=p.parseByte(field,i)-'0';require(digit>=0 && digit<=9,"Invalid decimal digit");value=value*10+digit;}
		require(value<=SCALE,"Fraction/moment exceeds one");return value;
	}
	private static long rendered(double x,ByteBuilder b){
		require(Double.isFinite(x) && x>=0 && x<=1,"Invalid moment");
		b.clear().append(x,12);
		if(b.length==1){
			require(b.array[0]=='0' || b.array[0]=='1',"Unexpected integral producer formatting");return (b.array[0]-'0')*SCALE;
		}
		require(b.length==14 && b.array[1]=='.',"Unexpected producer formatting");
		long value=b.array[0]-'0';for(int i=2;i<b.length;i++){value=value*10+b.array[i]-'0';}return value;
	}
	private static boolean sameRendered(double x,long expected,int forms,int field,ByteBuilder b){
		return rendered(x,b)==expected && (b.length==1)==((forms&(1<<field))!=0);
	}
	private static boolean equal(LineParser1 a,int ai,LineParser1 b,int bi){
		if(a.length(ai)!=b.length(bi)){return false;}for(int i=0;i<a.length(ai);i++){if(a.parseByte(ai,i)!=b.parseByte(bi,i)){return false;}}return true;
	}
	private static boolean matches(LineParser1 p,int field,ByteBuilder b){
		if(p.length(field)!=b.length){return false;}for(int i=0;i<b.length;i++){if(p.parseByte(field,i)!=b.array[i]){return false;}}return true;
	}
	private static long nonnegative(LineParser1 p,int field){final long n=p.parseLong(field);require(n>=0,"Negative count");return n;}
	private static long positive(String text){final long n=Long.parseLong(text);require(n>0,"Expected positive count");return n;}
	private static String take(HashMap<String,String> options,String name){final String value=options.remove(name);require(value!=null,"Missing argument: "+name);return value;}
	private static void header(ByteFile input,String header){require(Arrays.equals(input.nextLine(),header.getBytes(StandardCharsets.US_ASCII)),"Unexpected input header");}
	private static void require(boolean valid,String why){if(!valid){throw new IllegalArgumentException(why);}}
	private static double seconds(long ns){return ns*1e-9;}
	private static long cpuTime(){
		final java.lang.management.OperatingSystemMXBean bean=java.lang.management.ManagementFactory.getOperatingSystemMXBean();
		require(bean instanceof com.sun.management.OperatingSystemMXBean,"Native JVM process CPU timer required");
		final long time=((com.sun.management.OperatingSystemMXBean)bean).getProcessCpuTime();require(time>=0,"Process CPU timer unavailable");return time;
	}
	private void checkInputs()throws Exception{for(Snapshot s:inputs){s.check();}}
	private ByteStreamWriter writer(String name){return ByteStreamWriter.makeBSW(out.resolve(name).toString(),false,false,false);}
	private void writeText(String name,String text){final ByteStreamWriter w=writer(name);w.print(text);require(!w.poisonAndWait(),"Write failed: "+name);}
	private static final class Snapshot{
		Snapshot(Path path_)throws Exception{
			path=path_.toAbsolutePath().normalize();before=Files.readAttributes(path,BasicFileAttributes.class,LinkOption.NOFOLLOW_LINKS);
			require(before.isRegularFile() && Files.isReadable(path),"Expected readable ordinary input file");
		}
		void check()throws Exception{
			final BasicFileAttributes now=Files.readAttributes(path,BasicFileAttributes.class,LinkOption.NOFOLLOW_LINKS);
			require(now.isRegularFile() && now.size()==before.size() && now.lastModifiedTime().equals(before.lastModifiedTime()) &&
				String.valueOf(now.fileKey()).equals(String.valueOf(before.fileKey())),"Input changed: "+path);
		}
		final Path path;final BasicFileAttributes before;
	}
	private static final class Family{
		Family(int id_,int rank_,String label_,long expected_,long canonical_,long gc_,long ambiguous_,long ambiguousGenes_){
			id=id_;rank=rank_;label=label_;expected=expected_;expectedCanonical=canonical_;expectedGc=gc_;
			expectedAmbiguous=ambiguous_;expectedAmbiguousGenes=ambiguousGenes_;
		}
		final int id,rank;final String label;final long expected,expectedCanonical,expectedGc,expectedAmbiguous,expectedAmbiguousGenes;
		long expectedMean,expectedVariance=-1,expectedSd=-1,expectedMin,expectedMax;int integerForms;
		long n,visited,canonical,gc,ambiguous,ambiguousGenes;
		double mean,m2,sd,min=Double.POSITIVE_INFINITY,max=Double.NEGATIVE_INFINITY;
		final long[] excluded=new long[MULTIPLES.length];
	}
	private final String queryPath,familyPath,memberPath;private final Path out;
	private final long expectedQueries,expectedAssigned;private final int expectedFamilies;
	private final ArrayList<Family> families=new ArrayList<Family>();
	private final IntHashMap byId=new IntHashMap();private final ArrayList<Snapshot> inputs=new ArrayList<Snapshot>();
	private long unassignedUnknown;
	private static final long SCALE=1000000000000L;
	private static final String GRID="0,0.5,1,1.5,2,2.5,3,4,5,6";
	private static final double[] MULTIPLES={0,.5,1,1.5,2,2.5,3,4,5,6};
	private static final String[] SIZES={"all","1","2-4","5-9","10-49","50-99",">=100"};
	private static final String MEMBER_HEADER="input_index\tinput_id\tfamily_id\tstatus";
	private static final String QUERY_HEADER="input_index\tinput_id\tstable_family_id\tstatus\tgene_gc\ttail_gc\tcontig_gc\tgenome_gc\tgene_canonical\tgene_ambiguous\tcontig_canonical\tcontig_ambiguous\tgenome_canonical\tgenome_ambiguous";
	private static final String FAMILY_HEADER="stable_family_id\troster_rank\tfamily_id\texpected\tassigned\tmapped\tusable\tmissing_dna\tmissing_contig\tunusable\tambiguous_genes\tA\tC\tG\tT\tambiguous_bases\tmean_gc\tsample_variance\tsd\tmin_gc\tmax_gc";
}
