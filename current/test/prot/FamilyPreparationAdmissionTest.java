package test.prot;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

import ddl.DDLPentamerBlacklist;
import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import parse.LineParser1;
import prot.FamilyAlignmentMetric;
import prot.FamilyBlacklistBinding;
import prot.FamilyCalibration;
import prot.FamilyHeuristicMetadata;
import prot.FamilyMemberPreparation;
import prot.FamilyMemberPreparation.Prepared;
import prot.FamilyNegativeUpperCalibration;
import prot.FamilyNegativeUpperCalibration.MetricResult;
import prot.FamilyNegativeUpperCalibration.MetricWorker;
import prot.FamilyNegativeUpperCalibration.NegativePanel;
import prot.FamilyNegativeUpperCalibration.PackagedCalibration;
import prot.FamilyReferenceCodec;
import prot.FamilyReferenceExport;
import prot.FamilyTailCascade;
import prot.MemberLengthRange;
import prot.ProteinSubstitutionMatrix;
import shared.Shared;
import structures.ByteBuilder;

/** Exact admission-only API checks and one fixed throughput-only saved-data slice.
 * @author Sayu
 */
public final class FamilyPreparationAdmissionTest{

	public static void main(String[] args)throws Exception{
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg:args){final int equals=arg.indexOf('=');check(equals>0 && equals<arg.length()-1,"Use nonempty flag=value test arguments");
			check(options.put(arg.substring(0,equals),arg.substring(equals+1))==null,"Duplicate test argument");}
		check(options.containsKey("out") && options.containsKey("mode") && options.size()==(options.containsKey("in") ? 3 : 2),"Require mode=fixture|slice out=fresh_directory [in=accepted_recipe]");
		final Path out=Paths.get(options.get("out"));check(!Files.exists(out),"Preserve earlier fixture/measurement outputs");Files.createDirectory(out);
		Shared.AMINO_IN=true;Shared.setThreads(8);ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		if("fixture".equals(options.get("mode"))){fixture(out);}
		else{check("slice".equals(options.get("mode")) && options.containsKey("in"),"Only fixture or the single fixed slice is supported");slice(Paths.get(options.get("in")),out);}
	}

	private static void fixture(Path out)throws Exception{
		final String a=repeat('A',40),r=repeat('R',40),shortA=repeat('A',20);
		final Path negatives=out.resolve("negative.fa"),medians=out.resolve("medians.tsv"),hydro=out.resolve("hydro.tsv");
		write(negatives,">identical\n"+a+"\n>prefix\nR"+a+"\n>short\n"+shortA+"\n>other\n"+r+"\n>tie\n"+a+"\n");
		write(medians,"family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\nnormal\t1\t40\t40\t40\t1\nshortref\t1\t40\t40\t40\t1\n");
		write(hydro,"family_rank\tfamily_id\tN\thydrophobicity_mean\n0\tnormal\t1\t1.8\n1\tshortref\t1\t1.8\n");
		final FamilyAlignmentMetric metric=new FamilyAlignmentMetric(FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY,28);
		final FamilyMemberPreparation preparer=new FamilyMemberPreparation(ProteinKmerGeometry.exact(5),4,2,null,
			FamilyCalibration.MemberRetention.METRIC_SPECIFIC,true,new MemberLengthRange(40,99),ProteinPositionMode.FIXED,metric);
		final Prepared normal=preparer.prepare(0,"normal",new byte[][]{bytes(a)},bytes(a));
		final Prepared shortref=preparer.prepare(1,"shortref",new byte[][]{bytes(a)},bytes(shortA));
		final Prepared[] families={normal,shortref};
		final FamilyHeuristicMetadata metadata=new FamilyHeuristicMetadata(families,medians,pin(medians),hydro,pin(hydro));
		final FamilyTailCascade.ScalarBounds bounds=new FamilyTailCascade.ScalarBounds(new double[]{-4.5,-1,0,0,0},new double[]{4.5,1,1,1,20},1);
		final NegativePanel panel=FamilyNegativeUpperCalibration.loadPreparationPanel(negatives,5,20,41,bounds,ProteinKmerGeometry.exact(5),ProteinPositionMode.PLAIN,1,0,null);
		for(int floor:new int[]{1,1024}){
			final PackagedCalibration population=new PackagedCalibration(families,panel,metric,metadata,new int[]{0,1},new int[]{floor,floor});
			for(int top:new int[]{1,5}){
				final MetricWorker worker=population.newMetricWorker(top);
				for(int f=0;f<families.length;f++){
					final Trace old=new Trace(top),fast=new Trace(top);
					final MetricResult before=worker.scoreFamily(f,old::observe),after=worker.scoreFamilyForPreparation(f,fast::observe);
					same(before,after,old,fast,families[f],panel.rawCosineMaximum(families[f]));
				}
			}
		}
		// Independent witness: keep A spectra/consensus but deliberately bind the R sketch.
		// The identical-A query then has cosine1 yet no DDL hit; R is admitted with lower cosine.
		final Prepared rProfile=preparer.prepare(2,"r",new byte[][]{bytes(r)},bytes(r));
		System.arraycopy(rProfile.plain,0,normal.plain,0,normal.plain.length);
		final PackagedCalibration witness=new PackagedCalibration(new Prepared[]{normal},panel,metric,metadata,new int[]{0},new int[]{1});
		final MetricWorker worker=witness.newMetricWorker(5);final double[] admittedMax={0};
		final Trace old=new Trace(5),fast=new Trace(5);
		final MetricResult before=worker.scoreFamily(0,(f,q,h,c,o)->{old.observe(f,q,h,c,o);admittedMax[0]=Math.max(admittedMax[0],c);});
		final MetricResult after=worker.scoreFamilyForPreparation(0,fast::observe);
		check(Math.abs(before.cosineMaximum()-1)<1e-12 && before.cosineMaximum()>admittedMax[0] && before.admitted>0,"A rejected high-cosine query must remain in the public full maximum");
		same(before,after,old,fast,normal,panel.rawCosineMaximum(normal));
		boolean failed=false;
		try{worker.scoreFamilyForPreparation(0,(f,q,h,c,o)->{throw new IllegalArgumentException("observer abort");});}catch(IllegalArgumentException expected){failed=true;}
		check(failed && bits(worker.scoreFamily(0,null).cosineMaximum())==bits(before.cosineMaximum()),"Observer failure must propagate without poisoning the full-maximum worker");
		write(out.resolve("PASS"),"PREPARATION_ADMISSION_FIXTURE_PASS public_full_max=true rejected_max_witness=true selected_pairs=true states=true scientific_export=true raw400=true reset=true\n");
		System.out.print(new String(Files.readAllBytes(out.resolve("PASS")),StandardCharsets.US_ASCII));
	}

	private static void same(MetricResult old,MetricResult fast,Trace oldTrace,Trace fastTrace,Prepared family,double rawMaximum){
		check(old.hasCosineMaximum() && !fast.hasCosineMaximum() && Double.isFinite(fast.cosineMax) && fast.cosineMax==0,"The preparation-only maximum must be explicitly unavailable and finite");
		boolean refused=false;try{fast.cosineMaximum();}catch(IllegalStateException expected){refused=true;}check(refused,"Unmeasured combined maximum must refuse numeric access");
		check(old.scanned==fast.scanned && old.admitted==fast.admitted && old.selected==fast.selected && old.eligible==fast.eligible && old.familyRank==fast.familyRank && old.familyId.equals(fast.familyId),"Scientific count/identity parity must survive admission-first scoring");
		for(FamilyAlignmentMetric.State state:FamilyAlignmentMetric.State.values()){check(old.count(state)==fast.count(state),"Every selected alignment state must be conserved");}
		if(old.measured()){check(bits(old.identityUpper())==bits(fast.identityUpper()) && bits(old.blosumUpper())==bits(fast.blosumUpper()),"Measured identity/BLOSUM maxima must be bit-identical");}
		check(oldTrace.size==fastTrace.size && Arrays.equals(oldTrace.data,fastTrace.data),"Every selected ordinal, hit count, cosine and alignment observation must be identical");
		final ByteBuilder a=new ByteBuilder(),b=new ByteBuilder();
		final String schema=family.alignmentMetric.matrix.same(ProteinSubstitutionMatrix.BLOSUM62) ? FamilyReferenceCodec.ALIGNMENT_SCHEMA : FamilyReferenceCodec.MATRIX_SCHEMA;
		final FamilyReferenceCodec.Header header=new FamilyReferenceCodec.Header(schema,ProteinKmerGeometry.fromName(family.geometry).k(),family.width,0,1,1,0,FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING,0,family.alignmentMetric);
		FamilyReferenceCodec.appendFiveLines(FamilyReferenceExport.fromPrepared(family,rawMaximum,new MetricResult[]{old},0,0),a,header);
		FamilyReferenceCodec.appendFiveLines(FamilyReferenceExport.fromPrepared(family,rawMaximum,new MetricResult[]{fast},0,0),b,header);
		check(a.toString().equals(b.toString()),"Scientific exported reference bytes must not depend on the unused diagnostic");
	}

	/** First128 live families from original bucket0, all original k5 queries, eight workers.
	 * This fixed input is throughput-only, never a calibration or accuracy panel.
	 */
	@SuppressWarnings("unchecked")
	private static void slice(Path manifest,Path out)throws Exception{
		final FamilyCalibration.VerificationMode verification=FamilyCalibration.VerificationMode.STRUCTURAL;
		final Class<?> type=Class.forName("prot.FamilyPreparation");
		final Constructor<?> constructor=type.getDeclaredConstructor(Path.class,FamilyCalibration.VerificationMode.class);constructor.setAccessible(true);
		final Object recipe=constructor.newInstance(manifest,verification);
		final HashMap<String,String> o=(HashMap<String,String>)field(recipe,"options");
		final ProteinKmerGeometry geometry=(ProteinKmerGeometry)field(recipe,"geometry");
		final FamilyAlignmentMetric metric=(FamilyAlignmentMetric)field(recipe,"alignmentMetric");
		final MemberLengthRange range=(MemberLengthRange)field(recipe,"memberRange");
		final DDLPentamerBlacklist blacklist=FamilyBlacklistBinding.load(o,verification);
		final FamilyCalibration.MemberInputs source=new FamilyCalibration.MemberInputs(Paths.get(o.get("families")),o.get("familiessha80"),Paths.get(o.get("chunks")),o.get("chunkssha80"),Paths.get(o.get("consensus")),o.get("consensussha80"),Integer.parseInt(o.get("expectedfamilies")),Long.parseLong(o.get("expectedmembers")),Integer.parseInt(o.get("buckets")),Integer.parseInt(o.get("batches")),Integer.parseInt(o.get("maxreferenceaa")),range,verification);
		final FamilyCalibration.MemberInput[] bucket=source.loadBucket(0,Long.parseLong(o.get("maxbucketmembers")),Long.parseLong(o.get("maxbucketbytes")));
		Arrays.sort(bucket,(a,b)->Integer.compare(a.rank,b.rank));
		final Prepared[] prepared=new Prepared[128];final int[] ordinals=new int[128],floors=new int[128];
		final FamilyMemberPreparation preparer=new FamilyMemberPreparation(geometry,4,2,blacklist,FamilyCalibration.MemberRetention.METRIC_SPECIFIC,true,range,ProteinPositionMode.FIXED,metric,Integer.parseInt(o.get("width")));
		int n=0;for(FamilyCalibration.MemberInput row:bucket){if(row.members.length==0){continue;}prepared[n]=preparer.prepare(row.rank,row.id,row.members,row.consensus);ordinals[n]=row.rank;floors[n]=Integer.parseInt(o.get("ddlreject"));if(++n==128){break;}}
		check(n==128,"The fixed bucket must contain128 live families");
		final int count=Integer.parseInt(o.get("expectedfamilies"));final String[] ids=new String[count];final int[] ranks=new int[count];final long[] members=new long[count];
		final ByteFile roster=ByteFile.makeByteFile1(o.get("families"),false);final LineParser1 lp=new LineParser1('\t');int rows=0;
		try{roster.nextLine();for(byte[] line=roster.nextLine();line!=null;line=roster.nextLine()){lp.set(line);check(lp.terms()==3 && lp.parseInt(0)==rows && lp.parseLong(2)>0,"This fixed original cohort has a dense live roster");ids[rows]=lp.parseString(1);ranks[rows]=rows;members[rows++]=lp.parseLong(2);}}finally{check(!roster.close(),"Roster read must finish cleanly");}
		check(rows==count,"The original metadata roster must remain complete");
		final FamilyHeuristicMetadata metadata=new FamilyHeuristicMetadata(ids,ranks,members,Paths.get(o.get("lengthmedians")),o.get("lengthmedianssha80"),Paths.get(o.get("familyhydro")),o.get("familyhydrosha80"),verification);
		final ByteFile panels=ByteFile.makeByteFile1(o.get("negativepanels"),false);Path negative=null;int queries=0;
		try{panels.nextLine();for(byte[] line=panels.nextLine();line!=null;line=panels.nextLine()){lp.set(line);if(lp.termEquals("k5",0)){negative=Paths.get(lp.parseString(1));queries=lp.parseInt(3);}}}finally{check(!panels.close(),"Panel table read must finish cleanly");}
		check(negative!=null && queries==1246631,"The fixed measurement uses the original complete k5 panel");
		// The production helper is package-private; this fixture already reflects its recipe inputs.
		final java.lang.reflect.Method boundsMethod=FamilyTailCascade.class.getDeclaredMethod(
			"combinedBounds",HashMap.class,FamilyCalibration.VerificationMode.class);
		boundsMethod.setAccessible(true);
		final FamilyTailCascade.ScalarBounds bounds=(FamilyTailCascade.ScalarBounds)boundsMethod.invoke(null,o,verification);
		final NegativePanel panel=FamilyNegativeUpperCalibration.loadPreparationPanel(negative,queries,40,99,bounds,geometry,ProteinPositionMode.FIXED,4,2,blacklist,Integer.parseInt(o.get("width")));
		final PackagedCalibration population=new PackagedCalibration(prepared,panel,metric,metadata,ordinals,floors);
		final int top=Integer.parseInt(o.get("uppertop"));final MetricWorker warm=population.newMetricWorker(top);
		warm.scoreFamily(0,null);warm.scoreFamilyForPreparation(0,null);
		final Trace[][] traces={new Trace[128],new Trace[128]};final MetricResult[][] results={new MetricResult[128],new MetricResult[128]};final double[][] raw={new double[128],new double[128]};
		final java.lang.management.ThreadMXBean cpu=java.lang.management.ManagementFactory.getThreadMXBean();
		check(cpu.isCurrentThreadCpuTimeSupported() && cpu.isThreadCpuTimeEnabled(),"Fixed-slice worker CPU measurement must be supported and enabled");
		final long[][] workerCpu={new long[8],new long[8]};
		final ByteBuilder times=new ByteBuilder("mode\twall_seconds\tworker_cpu_seconds\tworker_selection_seconds\tworker_alignment_seconds\n");
		for(int mode=0;mode<2;mode++){
			final int selectedMode=mode;final AtomicInteger next=new AtomicInteger();final ExecutorService pool=Executors.newFixedThreadPool(8);final Future<?>[] tasks=new Future<?>[8];final long start=System.nanoTime();
			try{for(int w=0;w<8;w++){final int workerId=w;tasks[w]=pool.submit(()->{final long cpuStart=cpu.getCurrentThreadCpuTime();final MetricWorker worker=population.newMetricWorker(top);for(int f=next.getAndIncrement();f<128;f=next.getAndIncrement()){final Trace trace=new Trace(top);traces[selectedMode][f]=trace;results[selectedMode][f]=selectedMode==0 ? worker.scoreFamily(f,trace::observe) : worker.scoreFamilyForPreparation(f,trace::observe);raw[selectedMode][f]=panel.rawCosineMaximum(prepared[f]);}workerCpu[selectedMode][workerId]=cpu.getCurrentThreadCpuTime()-cpuStart;});}for(Future<?> task:tasks){task.get();}}finally{pool.shutdown();while(!pool.isTerminated()){pool.awaitTermination(1,java.util.concurrent.TimeUnit.SECONDS);}}
			long selection=0,alignment=0,usedCpu=0;for(MetricResult result:results[mode]){selection+=result.selectionNanos;alignment+=result.alignmentNanos;}for(long value:workerCpu[mode]){check(value>=0,"Worker CPU accounting must remain monotonic");usedCpu+=value;}
			times.append(mode==0 ? "full_max" : "preparation_only").tab().appendSlow((System.nanoTime()-start)*1e-9).tab().appendSlow(usedCpu*1e-9).tab().appendSlow(selection*1e-9).tab().appendSlow(alignment*1e-9).nl();
		}
		final ByteBuilder families=new ByteBuilder("rank\tid\tscanned\tadmitted\tselected\teligible\traw_cosine_maximum\n");
		for(int f=0;f<128;f++){check(bits(raw[0][f])==bits(raw[1][f]),"Required all-query raw400 maximum must remain bit-identical");same(results[0][f],results[1][f],traces[0][f],traces[1][f],prepared[f],raw[0][f]);final MetricResult result=results[0][f];families.append(result.familyRank).tab().append(result.familyId).tab().append(result.scanned).tab().append(result.admitted).tab().append(result.selected).tab().append(result.eligible).tab().appendSlow(raw[0][f]).nl();}
		write(out.resolve("timing.tsv"),times.toString());write(out.resolve("families.tsv"),families.toString());metadata.verify();
		write(out.resolve("PASS"),"PREPARATION_ADMISSION_SLICE_PASS families=128 panel=k5 queries=1246631 workers=8 selected_pairs=true scientific_export=true raw400=true calibration=false\n");
		System.out.print(times.toString());
	}

	private static Object field(Object value,String name)throws Exception{final Field field=value.getClass().getDeclaredField(name);field.setAccessible(true);return field.get(value);}
	private static final class Trace{
		Trace(int top){data=new long[top*8];}
		void observe(int family,int query,int hits,double cosine,FamilyAlignmentMetric.Observation value){data[size++]=query;data[size++]=hits;data[size++]=bits(cosine);data[size++]=value.state.ordinal();data[size++]=value.measured() ? 1 : 0;data[size++]=value.measured() ? bits(value.identityPct()) : 0;data[size++]=value.measured() ? bits(value.blosum()) : 0;data[size++]=family;}
		final long[] data;int size;
	}
	private static long bits(double value){return Double.doubleToLongBits(value);}
	private static String repeat(char acid,int length){final char[] text=new char[length];Arrays.fill(text,acid);return new String(text);}
	private static byte[] bytes(String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static String pin(Path file){return prot.DigestSuffix.file(file.toString());}
	private static void write(Path file,String value)throws Exception{Files.write(file,bytes(value));}
	private static void check(boolean valid,String reason){if(!valid){throw new AssertionError(reason);}}
}
