package test.prot;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import prot.*;
import prot.cluster.*;
import ddl.DDLFamilyProfiles;
import parse.LineParser1;
import structures.ByteBuilder;
/** Q45 software controls only: paired columns, package binding and caller compositions.
* The structural fixture builder is adapted from Qiqi's FamilyGcCallerTest.
* @author Sayu
*/
public final class FamilyPairedMinimumTest{
	/** Exercise literal gates, native construction and optional full caller composition. */
	public static void main(String[] args)throws Exception{
		check(args.length==1 || args.length==2 && args[1].equals("core"),"Fresh fixture root and optional core required");
		Path root=Paths.get(args[0]);
		Files.createDirectory(root);
		ProteinSubstitutionMatrix matrix=ProteinSubstitutionMatrix.BLOSUM62;
		byte[] w=new ProteinSequence("w",repeat('W',28)).enc;
		byte[] gapped=(repeat('m',22)+repeat('I',6)+repeat('D',6)).getBytes(StandardCharsets.US_ASCII);
		HbmPositionModel.TailResult gap=trace(w,w,gapped,matrix);
		TailAlignmentDecision old=gate(28,0,matrix).evaluate(w,w,gap),now=gate(28,28,matrix).evaluate(w,w,gap);
		check(old.accepted() && old.querySpan==28 && old.referenceSpan==28 && old.paired==22,"Old span-only gate admits gap-heavy28");
		check(now.reason==TailAlignmentDecision.Reason.PAIRED_SPAN_SHORT && now.paired==22 && Arrays.equals(now.path(),gapped),"Both-side gaps excluded; rejected trace retained");
		check(old.equals(new TailAlignmentGate(TailAlignmentGate.Mode.PIN_RIGHT_TRIM_LEFT,TailAlignmentGate.Denominator.ALIGNED_QUERY,100,28,0,0,null,matrix).evaluate(w,w,gap)),"Explicit zero matches old constructor exactly");
		System.out.println("GAP_HEAVY28 old="+old+" new="+now);
		byte[] a=new ProteinSequence("a",repeat('A',28)).enc,b=new ProteinSequence("s",repeat('S',28)).enc;
		TailAlignmentDecision mismatch=gate(28,28,matrix).evaluate(a,b,trace(a,b,repeat('m',28).getBytes(StandardCharsets.US_ASCII),matrix));
		check(mismatch.accepted() && mismatch.identities==0 && mismatch.paired==28,"Mismatches count as pairs");
		for(int n:new int[]{28,60,70}
		){
			byte[] seq=new ProteinSequence("edge",repeat('W',n)).enc;
			check(gate(n,n,matrix).decide(seq,seq).accepted(),"Inclusive gap-free paired boundary "+n);
			byte[] path=(repeat('m',n-1)+"ID").getBytes(StandardCharsets.US_ASCII);
			check(gate(n,n,matrix).evaluate(seq,seq,trace(seq,seq,path,matrix)).reason==TailAlignmentDecision.Reason.PAIRED_SPAN_SHORT,"Span n but pairs n-1 fails "+n);
			FamilyAlignmentMetric metric=new FamilyAlignmentMetric(FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY,n,matrix,n);
			FamilyAlignmentMetric.Observation out=new FamilyAlignmentMetric.Observation();
			ProteinSequence p=new ProteinSequence("p",repeat('W',n));
			metric.newWorker().measure(p,p,matrix.selfScore(p.enc),out);
			check(out.measured() && out.pairedColumns==n,"Measurement boundary "+n);
			for(String schema:new String[]{
				FamilyReferenceCodec.PAIRED_SCHEMA,FamilyReferenceCodec.PAIRED_MONOMER_SCHEMA
		}){
				FamilyReferenceCodec.Header h=header(schema,metric);
				byte[] bytes=FamilyReferenceCodec.appendHeader(h,new ByteBuilder()).toBytes();
				FamilyReferenceCodec.Header parsed=FamilyReferenceCodec.parseMetricBoundHeader(Arrays.copyOf(bytes,index(bytes,(byte)'\n')),new LineParser1('\t'));
				check(parsed.alignmentMetric.same(metric) && Arrays.equals(bytes,FamilyReferenceCodec.appendHeader(parsed,new ByteBuilder()).toBytes()),"Nonzero header roundtrip "+schema+" "+n);
			}
			rejects(()->header(FamilyReferenceCodec.MATRIX_SCHEMA,metric),"Old header cannot discard paired floor");
		}
		for(String schema:new String[]{
			FamilyReferenceCodec.ALIGNMENT_SCHEMA,FamilyReferenceCodec.MATRIX_SCHEMA,FamilyReferenceCodec.ALIGNMENT_MONOMER_SCHEMA,FamilyReferenceCodec.MATRIX_MONOMER_SCHEMA
		}){
			FamilyAlignmentMetric omitted=new FamilyAlignmentMetric(FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY,28,matrix);
			FamilyAlignmentMetric zero=new FamilyAlignmentMetric(omitted.mode,28,matrix,0);
			check(Arrays.equals(FamilyReferenceCodec.appendHeader(header(schema,omitted),new ByteBuilder()).toBytes(),FamilyReferenceCodec.appendHeader(header(schema,zero),new ByteBuilder()).toBytes()),"Off header bytes "+schema);
		}
		ProteinSequence query=new ProteinSequence("insert",repeat('W',14)+repeat('C',6)+repeat('W',14)),ref=new ProteinSequence("ref",repeat('W',28));
		FamilyAlignmentMetric.Observation obs=new FamilyAlignmentMetric.Observation();
		new FamilyAlignmentMetric(FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY,28,matrix,29).newWorker().measure(query,ref,matrix.selfScore(query.enc),obs);
		check(obs.state==FamilyAlignmentMetric.State.SPAN_TOO_SHORT && obs.pairedColumns<29,"Native metric rejects insufficient pairs without MEASURED");
		Path fasta=root.resolve("tiny.fa");
		write(fasta,">q\n"+repeat('W',100)+"\n");
		for(int n:new int[]{28,70}
		){
			DdlClust constructor=new DdlClust(new String[]{
				"in="+fasta,"out="+root.resolve("unused"+n+"#.faa"),"minlen=100","minspan="+n,
				"minpaired="+n,"matrix=blosum62","t=1"});
			java.lang.reflect.Field f=DdlClust.class.getDeclaredField("gate");
			f.setAccessible(true);
			check(((TailAlignmentGate)f.get(constructor)).minPaired==n,"Construction CLI paired floor "+n);
		}
		cake(root.resolve("cake"),fasta);
		if(args.length==2){
			write(root.resolve("PASS"),"FAMILY_PAIRED_CORE_TEST_PASS checks="+checks+"\n");
			System.out.println("FAMILY_PAIRED_CORE_TEST_PASS checks="+checks);
			return;
		}
		byte[] raw=sequence();
		FamilyTailCascade.LoadedClassifier off=fixture(root.resolve("off"),raw,new int[]{6,6,6},50,0,matrix);
		FamilyTailCascade.LoadedClassifier on=fixture(root.resolve("on"),raw,new int[]{6,6,6},50,65,matrix);
		check(off.newWorker().assign("off",raw).assigned(),"Off caller accepts exact fixture protein");
		check(!on.newWorker().assign("on",raw).assigned() && on.alignmentMinimumPaired()==65,"Caller rejects exact64 when paired floor65");
		FamilyScoringPackage recipe=new FamilyScoringPackage(root.resolve("on/classifier.config"),FamilyCalibration.VerificationMode.STRUCTURAL);
		check(recipe.load(1).alignmentMinimumPaired()==65,"Classifier recipe retains paired floor");
		HashMap<String,String> manifestOptions=recipe.resolvedOptions();
		manifestOptions.remove("minmemberaa");
		manifestOptions.remove("maxmemberaa");
		for(String[] pair:new String[][]{{"heuristic","heuristicarm"},
			{"lengthmedians","familymedians"},{"lengthmedianssha80","familymedianssha80"}}){
			manifestOptions.put(pair[1],manifestOptions.remove(pair[0]));
		}
		manifestOptions.put("minqueryaa","100");
		manifestOptions.put("maxqueryaa","100");
		manifestOptions.put("marginscale","1");
		StringBuilder manifest=new StringBuilder("#bbgap_tail_cascade\t2\n");
		for(String key:new TreeSet<String>(manifestOptions.keySet())){
			manifest.append(key).append('\t').append(manifestOptions.get(key)).append('\n');
		}
		Path manifestFile=root.resolve("caller.tsv");
		write(manifestFile,manifest.toString());
		check(prok.StopTailCascadeResources.readFixed100Structural(manifestFile.toString(),"unmeasured").load(1).alignmentMinimumPaired()==65,"Actual caller manifest propagates paired floor");
		HashMap<String,String> opts=recipe.resolvedOptions();
		opts.remove("heuristic");
		opts.remove("lengthmedians");
		opts.remove("lengthmedianssha80");
		opts.remove("familyhydro");
		opts.remove("familyhydrosha80");
		opts.remove("minqueryaa");
		opts.remove("maxqueryaa");
		opts.put("alignmentminpaired","28");
		rejects(()->FamilyTailCascade.loadCallableStructural(opts,1,new ddl.DDLFamilySearch.StructuralInputs()),"Callable/header mismatch rejected");
		opts.remove("alignmentminpaired");
		rejects(()->FamilyTailCascade.loadCallableStructural(opts,1,new ddl.DDLFamilySearch.StructuralInputs()),"Omitted floor cannot load paired header");
		FamilyTailCascade.LoadedClassifier length=on.withLengthPolicy(policy());
		check(length.alignmentMinimumPaired()==65 && length.alignmentMinimumSpan()==28 && !length.newWorker().assign("length",raw).assigned(),"Length-policy composition retains paired floor");
		int[][] scores=new int[20][20];
		for(int i=0;i<20;i++){
			for(int j=0;j<20;j++){scores[i][j]=4*Blosum62.score((byte)i,(byte)j);}
		}
		ProteinSubstitutionMatrix pairwise=ProteinSubstitutionMatrix.fromScores("pairwise40",4,16,scores),other=ProteinSubstitutionMatrix.fromScores("pairwise35",4,16,scores);
		FamilyTailCascade.LoadedClassifier swapped=fixture(root.resolve("swap"),raw,new int[]{6,6,6},50,65,pairwise).withUncalibratedMatrix(other);
		check(swapped.alignmentMinimumPaired()==65 && !swapped.newWorker().assign("swapped",raw).assigned(),"Matrix replacement retains paired floor");
		write(root.resolve("PASS"),"FAMILY_PAIRED_MINIMUM_TEST_PASS checks="+checks+"\n");
		System.out.println("FAMILY_PAIRED_MINIMUM_TEST_PASS checks="+checks);
	}
	/** Run a real resident structural layer, retaining its replayed paired option. */
	private static void cake(Path result,Path fasta)throws Exception{
		ProteinSequence sequence=new ProteinSequence("seed",repeat('W',100));
		ddl.FamilySketchBuilder builder=new ddl.FamilySketchBuilder(64,2,7,false,0,0,null);
		Family family=new Family(0,"seed",0,sequence.enc);
		family.refreshStats(new ProteinSequence[]{sequence},builder);
		Path codec=result.resolveSibling("tiny.pddl");
		ClusterSketchCodec.write(codec.toString(),false,new Family[]{family},2,64,7,false,0,0);
		DdlClust.main(new String[]{
			"mode=cake","memberranking=f","in="+fasta,"codecs="+codec,"outdir="+result,
			"familybudget=1","initialfraction=1","passes=1","fractions=1","lambda=0","t=1",
			"k=2","width=64","sketchseed=7","positional=f","partitions=0","overlap=0",
			"minlen=100","minspan=28","minpaired=28","matrix=blosum62","maxalignlen=100","minid=1","minblosumratio=0",
			"mergeminid=1","mergeminratio=0","trimleft=f","consensustrimleft=f","quality=null",
			"passhits=0","greedyhits=0","maxhydrodelta=9","greedymaxhydrodelta=9","fissiontopfraction=0","fusionbottomfraction=0"
		});
		check(Files.readAllLines(result.resolve("layer1/run.config"),StandardCharsets.US_ASCII).contains("minpaired=28"),"Resident structural replay retains paired28");
		check(Files.isRegularFile(result.resolve("final/families.pddl.gz")),"Real structural layer completes final model output");
	}
	/** Java8-safe repeated fixture residues and trace operators. */
	private static String repeat(char c,int n){
		char[] a=new char[n];
		Arrays.fill(a,c);
		return new String(a);
	}
	/** Synthetic policy certificate used only to exercise metric composition. */
	static FamilyLengthPolicy policy(){
		FamilyLengthPolicy.Requirement[] rows=new FamilyLengthPolicy.Requirement[8];
		for(int i=0;i<8;i++){
			rows[i]=new FamilyLengthPolicy.Requirement(FamilyLengthPolicy.view(i),FamilyLengthPolicy.cohort(i),FamilyLengthPolicy.expectedQueries(i),FamilyLengthPolicy.expectedBaseline(i),Math.toIntExact(FamilyLengthFit.allowed100ppm(FamilyLengthPolicy.expectedQueries(i))),0,0,0,0,"/fixture/input","00000000000000000000","/fixture/report");
		}
		return new FamilyLengthPolicy("/fixture/reference","unmeasured","0000000000000000000000000000000000000000",rows);
	}
	/** Build a minimal explicit codec header. */
	static FamilyReferenceCodec.Header header(String schema,FamilyAlignmentMetric metric){
		return new FamilyReferenceCodec.Header(schema,5,1024,0,0,0,0,FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING,0,metric);
	}
	/** Use permissive identity/score gates so only span and pair floors decide. */
	static TailAlignmentGate gate(int span,int paired,ProteinSubstitutionMatrix matrix){
		return new TailAlignmentGate(TailAlignmentGate.Mode.PIN_RIGHT_TRIM_LEFT,TailAlignmentGate.Denominator.ALIGNED_QUERY,100,span,0,0,null,matrix,paired);
	}
	/** Score a literal recorded trace independently for gate validation. */
	static HbmPositionModel.TailResult trace(byte[] q,byte[] r,byte[] path,ProteinSubstitutionMatrix matrix){
		int score=0,i=0,j=0;
		for(byte op:path){
			if(op=='m'){
				score+=matrix.score(q[i++],r[j++]);
			}
			else{
				score-=matrix.linearGap;
				if(op=='I'){i++;}
				else{j++;}
			}
		}
		return new HbmPositionModel.TailResult(score*HbmPositionModel.SCALE,0,q.length-1,0,r.length-1,path);
	}
	/** Locate the first serialized header line. */
	static int index(byte[] a,byte b){
		for(int i=0;i<a.length;i++){
			if(a[i]==b){return i;}
		}
		throw new AssertionError();
	}
	interface Throwing{
		void run()throws Exception;
	}
	/** Require an explicit invalid-binding failure. */
	static void rejects(Throwing run,String why)throws Exception{
		try{
			run.run();
		}
		catch(IllegalArgumentException expected){
			check(true,why);
			return;
		}
		throw new AssertionError(why);
	}
	/** Qiqi-derived structural classifier with explicit matrix and paired binding. */
	public static FamilyTailCascade.LoadedClassifier fixture(final Path root,final byte[] raw,final int[] hits,final int top,final int paired,final ProteinSubstitutionMatrix matrix)throws Exception{
		check(hits!=null && hits.length>=3,"The fixture needs at least three families");
		for(int h:hits){check(h>=0,"Explicit nonnegative self-hit counts required");}
		Files.createDirectory(root);
		final char[] sketch=new char[1024]; DDLFamilyProfiles.sketchPositional(raw,5,0,4,2,sketch);
		final FamilyAlignmentMetric metric=new FamilyAlignmentMetric(FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY,28,matrix,paired);
		final FamilyReferenceCodec.Header header=new FamilyReferenceCodec.Header(paired>0 ? FamilyReferenceCodec.PAIRED_SCHEMA : matrix.same(ProteinSubstitutionMatrix.BLOSUM62) ? FamilyReferenceCodec.ALIGNMENT_SCHEMA : FamilyReferenceCodec.MATRIX_SCHEMA,
			5,1024,0,hits.length,hits.length,0,FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING,3.8906,metric);
		final ByteBuilder refs=FamilyReferenceCodec.appendHeader(header,new ByteBuilder());
		final ByteBuilder profiles=new ByteBuilder().append(DDLFamilyProfiles.HEADER16).nl();
		for(int f=0; f<hits.length; f++){
			final char[] row=sketch.clone(); int available=0;
			for(int bin=0; bin<row.length; bin++){if(row[bin]!=0 && ++available>hits[f]){row[bin]=0;}}
			check(available>=hits[f],"The fixture query must supply at least five distinct positional registers");
			final String id="f"+f; final float[] spectrum=new float[400]; spectrum[0]=1;
			FamilyReferenceCodec.appendFiveLines(new FamilyReference(id,f,30,1920,1800,new ProteinSequence(id,raw),row,spectrum,
				0,30,0,1.1f,99,1,raw,0,new double[10],FamilyReference.AlignmentCalibration.ENABLED),refs,header);
			profiles.append("tail\t").append(id).append("\t1024\t30\t1920\t1800\tnonempty\t");
			for(int field=0; field<3; field++){
				if(field>0){profiles.tab();}
				for(int bin=0; bin<row.length; bin++){
					if(bin>0){profiles.append(',');} profiles.append(field==0 ? (int)row[bin] : row[bin]==0 ? 0 : 1);
				}
			} profiles.nl();
		}
		final Path ref=root.resolve("references.tsv"),data=root.resolve("profiles.tsv"),summary=root.resolve("summary.tsv"),
			list=root.resolve("list.tsv"),bounds=root.resolve("bounds.tsv"),medians=root.resolve("medians.tsv"),hydro=root.resolve("hydro.tsv");
		Files.write(ref,refs.toBytes()); Files.write(data,profiles.toBytes());
		write(summary,"namespace\tk\tseed\texponent\tregisterbits\tsketch\tpartitions\toverlap\ttag_bits\tsegment_words\n"+
			"tail\t5\t0\t6\t16\t"+DDLFamilyProfiles.SKETCH_POSITIONAL+"\t4\t2\t"+
			DDLFamilyProfiles.positionalTagBits(4)+"\t"+DDLFamilyProfiles.positionalSegmentWords(4)+"\n");
		write(list,data+"\t"+unmeasured(data)+"\t"+unmeasured(summary)+"\n");
		final StringBuilder scalar=new StringBuilder("scalar\trows\tmin\tmax\tmin_family\tmax_family\n");
		for(String name:new String[]{"hydrophobicity","net_charge","aromaticity","fraction_small","pI_tail"}){
			scalar.append(name).append("\t79565\t-20\t20\tlo\thi\n");
		} write(bounds,scalar.toString());
		final StringBuilder ms=new StringBuilder("family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\n");
		final StringBuilder hs=new StringBuilder("family_rank\tfamily_id\tN\thydrophobicity_mean\n");
		for(int f=0; f<hits.length; f++){ms.append('f').append(f).append("\t30\t64\t64\t64\t30\n");hs.append(f).append("\tf").append(f).append("\t30\t1.8\n");}
		write(medians,ms.toString()); write(hydro,hs.toString());
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String entry:new String[]{"top=50","k=5","seed=0","width=1024","ddlreject=6","ddlaccept=1025",
			"positional=t","partitions=4","overlap=2","indexearlyreject=f","cosine=ranking","rankingfeatures=combined25",
			"scalargates=t","scalarz=3.8906","identityacceptfloor=30","blosumacceptfloor=0.145","scalarscale=0.1",
			"minmemberaa=40","maxmemberaa=99","alignmentmode=pin_right_trim_left_aligned_query","alignmentminspan=28"}){
			final int split=entry.indexOf('='); options.put(entry.substring(0,split),entry.substring(split+1));
		}
		options.put("matrix",matrix.name);options.put("matrixbinding",matrix.binding());
		if(paired>0){options.put("alignmentminpaired",Integer.toString(paired));}
		options.put("familyreferences",ref.toString()); options.put("familyreferencessha80",unmeasured(ref));
		options.put("top",Integer.toString(top));
		options.put("ddlprofiles",list.toString()); options.put("ddlprofilessha80",unmeasured(list));
		options.put("scalarbounds",bounds.toString()); options.put("scalarboundssha80",unmeasured(bounds));
		final FamilyTailCascade.LoadedClassifier base=FamilyTailCascade.loadCallableStructural(options,1,new ddl.DDLFamilySearch.StructuralInputs());
		final StringBuilder recipe=new StringBuilder();
		for(String key:new java.util.TreeSet<String>(options.keySet())){recipe.append(key).append('=').append(options.get(key)).append('\n');}
		recipe.append("heuristic=hitfourth_h020_c150\nminqueryaa=40\nmaxqueryaa=99\nlengthmedians=").append(medians)
			.append("\nlengthmedianssha80=").append(unmeasured(medians)).append("\nfamilyhydro=").append(hydro)
			.append("\nfamilyhydrosha80=").append(unmeasured(hydro)).append('\n');
		write(root.resolve("classifier.config"),recipe.toString());
		return base.withHeuristicRanking(new FamilyHeuristicScore(base,medians,unmeasured(medians),hydro,unmeasured(hydro),40,99,"hitfourth_h020_c150",prot.FamilyCalibration.VerificationMode.STRUCTURAL));
	}
	/** Deterministic small protein with enough positional registers for the index. */
	public static byte[] sequence(){
		final byte[] alphabet="ACDEFGHIKLMNPQRSTVWY".getBytes(StandardCharsets.US_ASCII),raw=new byte[64];
		for(int i=0; i<raw.length; i++){raw[i]=alphabet[(i*7+i/20)%alphabet.length];} return raw;
	}

	/** Synthetic structural inputs require a syntactic inherited label in scalar
	 * and heuristic loaders. This is a dummy, never a measured digest or provenance.
	 * Matches FamilyEvalStructuralTest's native no-hash fixture convention. */
	private static String unmeasured(Path path){return "00000000000000000000";}
	private static void write(Path path,String value)throws Exception{Files.write(path,value.getBytes(StandardCharsets.US_ASCII));}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}checks++;}
	private static int checks;
}
