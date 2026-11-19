package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import ddl.DDLFamilyProfiles;
import ddl.DDLFamilySearch;
import ddl.DDLPentamerBlacklist;
import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import template.ThreadWaiter;

/** Per-family negative upper-gate selection for family-tail cascade calibration.
 * Tracks every cosine control and each family's hardest query ordinals for forced DP.
 * @author Yelan, Keqing
 */
public final class FamilyNegativeUpperCalibration{

	private FamilyNegativeUpperCalibration(){}

	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args,FamilyNegativeUpperCalibration.class,false);
			try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	public static void run(String[] args) throws Exception{
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			if(pair.length!=2){
				throw new IllegalArgumentException("Require flag=value arguments");
			}
			final String key=pair[0].toLowerCase(Locale.ROOT);
			if(!OPTIONS.contains(key) || options.put(key,pair[1])!=null){
				throw new IllegalArgumentException("Unknown or duplicate option: "+key);
			}
		}
		final int top=integer(options,"top"), k=integer(options,"k"), width=integer(options,"width");
		if(options.containsKey("aabits") && options.containsKey("aapartition")){
			throw new IllegalArgumentException("Use either aabits or aapartition, not both");
		}
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(k, options.containsKey("aabits") ? options.get("aabits") : options.get("aapartition"));
		final long seed=number(options,"seed"), maxReads=number(options,"maxreads");
		final int threads=threads(options), shard=integer(options,"shard",0), shards=integer(options,"shards",1);
		if(top<1 || geometry.k()<4 || geometry.k()>6 || width!=1024 || seed!=0 || maxReads<1 || shard<0 || shard>=shards){
			throw new IllegalArgumentException("Upper calibration requires top>=1, k=4..6, width=1024, seed=0, maxreads>=1");
		}
		if(options.containsKey("familyreferences")){
			runPackaged(options,args,geometry,top,threads,shard,shards,Math.toIntExact(maxReads)); return;
		}
		if(options.containsKey("consensus")!=options.containsKey("consensussha80")){
			throw new IllegalArgumentException("consensus and consensussha80 are required together");
		}
		if(options.containsKey("queryids") || options.containsKey("queryidssha80")){
			throw new IllegalArgumentException("queryids and queryidssha80 require familyreferences");
		}
		if(containsAny(options,PACKAGED_REPLACEMENT_OPTIONS)){
			throw new IllegalArgumentException("DDL replacement and blacklist options require familyreferences");
		}
		final String referencesHash=required(options,"referencessha80"), modelsHash=required(options,"modelssha80");
		final String negativesHash=required(options,"negativessha80");
		final Path references=input(options,"references"), models=input(options,"models"), negatives=input(options,"negatives");
		final Path consensus=options.containsKey("consensus") ? input(options,"consensus") : null;
		final Path out=Paths.get(required(options,"out")).toAbsolutePath().normalize();
		if(Files.exists(out) || !Files.isDirectory(out.getParent())){
			throw new IllegalArgumentException("Output must be fresh under an existing directory");
		}
		Shared.AMINO_IN=true; Shared.setThreads(threads); ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final long started=System.nanoTime();
		final ReferenceFiles referenceFiles=loadReferenceManifest(references);
		final DDLFamilySearch.Index index=loadIndex(referenceFiles,geometry,seed,width,out.getParent());
		final Models fixed=consensus==null ? loadModels(models,index,integer(options,"maxreferenceaa")) :
			loadModels(models,consensus,index,integer(options,"maxreferenceaa"));
		final LinkedSpectra spectra=loadSpectra(referenceFiles.spectra,index.names);
		final ArrayList<Negative> queries=loadNegatives(negatives,Math.toIntExact(maxReads),
			integer(options,"minnegativeaa"),integer(options,"maxnegativeaa"),geometry,seed,width,null,false,1,0,null,false);
		final double[] norms=new double[spectra.values.length];
		for(int i=0;i<norms.length;i++){norms[i]=ProteinSpectrum.norm(spectra.values[i]);}
		final long loaded=System.nanoTime();
		final FamilyResult[] results=scoreFamilies(index,fixed.values,spectra.values,norms,queries,top,threads,shard,shards);
		long cosineCalls=0,alignments=0;
		for(FamilyResult result : results){
			if(result!=null){cosineCalls+=result.observations; alignments+=result.aligned;}
		}
		Files.createDirectory(out);
		final ByteStreamWriter upper=writer(out.resolve("uppers.tsv"));
		final ByteBuilder line=new ByteBuilder();
		try{
			upper.println("roster_rank\tfamily_id\tnegative_count\tcosine_upper\taligned_count\tidentity_pct_upper\tblosum_a_sq_upper");
			for(int family=0; family<index.families(); family++){
				final FamilyResult result=results[family];
				if(result==null){continue;}
				line.clear().append(fixed.values[family].rank).tab().append(fixed.values[family].id).tab()
					.append(result.observations).tab().append(Double.toString(result.cosineMax)).tab()
					.append(result.aligned).tab().append(Double.toString(result.identity)).tab()
					.append(Double.toString(result.blosum));
				upper.println(line);
			}
		}finally{if(upper.poisonAndWait()){throw new IllegalStateException("Upper output writer failed");}}
		writeTimings(out,results,fixed.values,index.families(),queries.size(),cosineCalls,alignments,top,threads,shard,shards,started,loaded);
		write(out.resolve("recipe.txt"),String.join("\n",args)+"\nschema=family_negative_upper_calibration_v1\n"+
			"score=ddl_matches_plus_dimer_cosine\nnegative_tie=ordinal_ascending\nthreads="+threads+
			"\nshard="+shard+"\nshards="+shards+"\n");
		pin(references,referencesHash); pin(models,modelsHash); pin(negatives,negativesHash);
		if(consensus!=null){pin(consensus,required(options,"consensussha80"));}
		loadReferenceManifest(references);
		write(out.resolve("PASS"),"FAMILY_NEGATIVE_UPPER_CALIBRATION_PASS\n");
		final StringBuilder pins=new StringBuilder("file\tsha80\n");
		for(String name : new String[]{"uppers.tsv","summary.tsv","recipe.txt","timing.tsv","PASS"}){
			pins.append(name).append('\t').append(DigestSuffix.file(out.resolve(name).toString())).append('\n');
		}
		write(out.resolve("result.sha80.tsv"),pins.toString());
		System.out.println("FAMILY_NEGATIVE_UPPER_CALIBRATION_PASS families="+index.families()+
			" negatives="+queries.size()+" alignments="+alignments);
	}
	/** Full family-side topN over every negative; no query-side shortlist or early gate. */
	private static void runPackaged(HashMap<String,String> o,String[] argv,ProteinKmerGeometry geometry,int top,int threads,int shard,int shards,int n) throws Exception{
		if(o.containsKey("consensus") || o.containsKey("consensussha80")){
			throw new IllegalArgumentException("consensus= applies only to legacy model-table upper calibration");
		}
		if(top>n || n>MAX_PANEL_QUERIES || o.containsKey("references") || o.containsKey("models")){
			throw new IllegalArgumentException("Packaged calibration requires top<=queries<="+MAX_PANEL_QUERIES+" and no legacy inputs");
		}
		if(o.containsKey("queryids")!=o.containsKey("queryidssha80")){
			throw new IllegalArgumentException("queryids and queryidssha80 are required together");
		}
		if(o.containsKey("queryprefix")==o.containsKey("queryids")){
			throw new IllegalArgumentException("Packaged calibration requires exactly one of queryprefix= or queryids=");
		}
		final Path refs=input(o,"familyreferences"),negatives=input(o,"negatives"),boundsFile=input(o,"scalarbounds");
		final FamilyTailCascade.ScalarBounds bounds=FamilyTailCascade.combinedBounds(o);
		final ProteinPositionMode positionMode=ProteinPositionMode.resolve(o.get("positionmode"),null,
			o.containsKey("positional") ? Boolean.valueOf(optionalBoolean(o,"positional",false)) : null,ProteinPositionMode.PLAIN);
		final boolean positional=positionMode.positional();
		final int partitions=o.containsKey("partitions") ? integer(o,"partitions") : 4;
		final int overlap=o.containsKey("overlap") ? integer(o,"overlap") : 8;
		final boolean replacement=positional || o.containsKey("ddlprofiles");
		if(!geometry.isExact() && !replacement){throw new IllegalArgumentException("Irregular amino-acid geometry requires replacement profiles");}
		if(!positional && (o.containsKey("partitions") || o.containsKey("overlap")) ||
				o.containsKey("ddlprofiles")!=o.containsKey("ddlprofilessha80")){
			throw new IllegalArgumentException("DDL replacement requires path+sha80; partitions/overlap require positional=t");
		}
		if(positional){DDLFamilyProfiles.validatePositionalLayout(partitions,overlap,geometry);}
		if(o.containsKey("blacklist")!=o.containsKey("blacklistsha80") ||
				(o.containsKey("blacklist") && !replacement)){
			throw new IllegalArgumentException("Experimental pentamer blacklist requires path+sha80 and matching profile replacement");
		}
		final Path ddlProfiles=replacement ? input(o,"ddlprofiles") : null;
		final DDLPentamerBlacklist blacklist=o.containsKey("blacklist") ?
			DDLPentamerBlacklist.load(input(o,"blacklist"),required(o,"blacklistsha80")) : null;
		final String prefix=o.get("queryprefix");
		if(prefix!=null && !prefix.matches("[A-Za-z0-9_]+")){throw new IllegalArgumentException("Unsafe query prefix");}
		final Path out=Paths.get(required(o,"out")).toAbsolutePath().normalize();
		if(Files.exists(out) || !Files.isDirectory(out.getParent())){throw new IllegalArgumentException("Output must be fresh");}
		Shared.AMINO_IN=true; Shared.setThreads(threads); ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final long start=System.nanoTime();
		final FamilyReferenceLoader.Loaded loaded=FamilyReferenceLoader.load(refs.toString(),threads,integer(o,"maxfamilies"));
		if(loaded.all.length!=integer(o,"maxfamilies")){throw new IllegalArgumentException("Exact parent family census required");}
		final char[][] fullReplacementRows=replacement ?
			FamilyTailCascade.loadReplacementRows(loaded.live,ddlProfiles,geometry,positionMode,partitions,overlap,blacklist) : null;
		final ArrayList<FamilyReference> selected=new ArrayList<FamilyReference>();
		final ArrayList<char[]> replacementRows=replacement ? new ArrayList<char[]>() : null;
		int empty=0;
		for(int i=0;i<loaded.live.length;i++){
			final FamilyReference f=loaded.live[i];
			if(f.rank%shards==shard){
				selected.add(f);
				if(replacementRows!=null){replacementRows.add(fullReplacementRows[i]);}
			}
		}
		for(FamilyReference f:loaded.all){if(f.rank%shards==shard && f.members==0){empty++;}}
		if(selected.isEmpty()){throw new IllegalArgumentException("Selected family shard must contain live families");}
		final FamilyReference[] families=selected.toArray(new FamilyReference[selected.size()]);
		final NegativePanel queries=loadPackagedPanel(negatives,n,integer(o,"minnegativeaa"),integer(o,"maxnegativeaa"),
			bounds,geometry,positionMode,partitions,overlap,blacklist);
		final PackagedCalibration calibration=new PackagedCalibration(families,
			replacementRows==null ? null : replacementRows.toArray(new char[replacementRows.size()][]),queries);
		final Path queryIds=prefix==null ? input(o,"queryids") : null;
		if(prefix==null){validateQueryIds(queryIds,queries);}
		else{for(int q=0;q<n;q++){if(!queries.queryId(q).equals(prefix+q)){throw new IllegalArgumentException("Negative query IDs must be canonical ordinals");}}}
		final long loadedAt=System.nanoTime(); Files.createDirectory(out);
		final AtomicInteger next=new AtomicInteger(); final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final FamilyResult[] results=new FamilyResult[families.length]; final ArrayList<Thread> workers=new ArrayList<Thread>();
		final int workerCount=Math.min(threads,families.length);
		for(int worker=0;worker<workerCount;worker++){
			final int w=worker;
			workers.add(new Thread(()->{
				ByteStreamWriter raw=null;
				try{
					raw=writer(out.resolve("raw_pairs_"+w+".tsv")); raw.println(PACKAGED_RAW_HEADER);
					final ByteStreamWriter output=raw; final ByteBuilder line=new ByteBuilder();
					final PackagedWorker scorer=calibration.newWorker(top);
					final PairObserver observer=(f,ordinal,ddl,cosine,identity,blosum)->{
						final FamilyReference family=families[f];
						line.clear().append(family.rank).tab().append(family.id).tab().append(ordinal).tab().append(queries.queryId(ordinal)).tab()
							.append(ddl).tab().append(Double.toString(cosine)).tab().append(Double.toString(identity)).tab().append(Double.toString(blosum)); output.println(line);
					};
					for(int f=next.getAndIncrement();f<families.length && failure.get()==null;f=next.getAndIncrement()){
						results[f]=scorer.scoreFamily(f,observer);
					}
				}catch(Throwable t){failure.compareAndSet(null,t);}
				finally{if(raw!=null && raw.poisonAndWait()){failure.compareAndSet(null,new IllegalStateException("Raw pair writer failed"));}}
			},"family-side-top500-"+worker));
		}
		startAndJoin(workers,failure); if(failure.get()!=null){throw new RuntimeException("Family-side calibration worker failed",failure.get());}
		final ByteBuilder upper=new ByteBuilder(FamilyReferenceExport.SUPPLEMENTAL_HEADER+"\n"),coverage=new ByteBuilder("roster_rank\tfamily_id\tnegative_count\taligned_count\tstatus\n");
		long pairs=0;
		for(int f=0;f<families.length;f++){
			final FamilyResult r=results[f]; final FamilyReference family=families[f];
			if(r==null || r.observations!=n || r.aligned!=top){throw new IllegalStateException("Incomplete family-side calibration");}
			pairs+=r.aligned; upper.append(family.rank).tab().append(family.id).tab().append(r.aligned).tab().append(Double.toString(r.identity)).tab().append(Double.toString(r.blosum)).nl();
			coverage.append(family.rank).tab().append(family.id).tab().append(n).tab().append(r.aligned).tab().append("MEASURED").nl();
		}
		for(FamilyReference f:loaded.all){if(f.rank%shards==shard && f.members==0){coverage.append(f.rank).tab().append(f.id).tab().append(0).tab().append(0).tab().append("NO_MEMBERS").nl();}}
		write(out.resolve("uppers.tsv"),upper.toString()); write(out.resolve("coverage.tsv"),coverage.toString());
		write(out.resolve("summary.tsv"),"families\tlive\tempty\tqueries\tcomparisons\talignments\ttop\tshard\tshards\tload_seconds\ttotal_seconds\n"+
			(families.length+empty)+"\t"+families.length+"\t"+empty+"\t"+n+"\t"+((long)n*families.length)+"\t"+pairs+"\t"+top+"\t"+shard+"\t"+shards+"\t"+
			((loadedAt-start)*1e-9)+"\t"+((System.nanoTime()-start)*1e-9)+"\n");
		write(out.resolve("recipe.txt"),String.join("\n",argv)+"\nselection=ALL_queries_per_live_family;DDL_plus_combined25_desc_then_query_ordinal_asc\n"+
			"rank_key=FamilyTailCascade.combinedRankingScore;shared_with_deployed_COMBINED_COMPARATOR\nscalar_autofail=false\nraw_scores=original_double\n"+
			"ddl_replacement="+replacement+"\npositional="+positional+"\npartitions="+(positional ? partitions : 1)+"\noverlap="+(positional ? overlap : 0)+"\n"+
			"blacklist="+(blacklist==null ? "none" : blacklist.sha80()+";words="+blacklist.size())+"\n"+
			"empty_family=NO_MEMBERS;no_consensus_no_alignment\n");
		pin(refs,required(o,"familyreferencessha80")); pin(negatives,required(o,"negativessha80")); pin(boundsFile,required(o,"scalarboundssha80"));
		if(queryIds!=null){pin(queryIds,required(o,"queryidssha80"));}
		write(out.resolve("PASS"),"FAMILY_NEGATIVE_PACKAGED_PASS\n");
		final StringBuilder manifest=new StringBuilder("file\tsha80\n");
		for(String name:new String[]{"uppers.tsv","coverage.tsv","summary.tsv","recipe.txt","PASS"}){manifest.append(name).append('\t').append(DigestSuffix.file(out.resolve(name).toString())).append('\n');}
		for(int w=0;w<workerCount;w++){final String name="raw_pairs_"+w+".tsv"; manifest.append(name).append('\t').append(DigestSuffix.file(out.resolve(name).toString())).append('\n');}
		write(out.resolve("result.sha80.tsv"),manifest.toString()); System.out.println("FAMILY_NEGATIVE_PACKAGED_PASS live="+families.length+" empty="+empty+" pairs="+pairs);
	}

	/** Loads owned, immutable query features in FASTA order using the existing numerical path.
	 * The caller configures amino-acid input and validates file pins and query identities.
	 * This panel is for packaged combined25 selection, not the legacy raw-dimer cosine gate.
	 */
	public static NegativePanel loadPackagedPanel(Path file, int expected, int minLength, int maxLength,
			FamilyTailCascade.ScalarBounds bounds, boolean positional, int partitions, int overlap, DDLPentamerBlacklist blacklist){
		return loadPackagedPanel(file, expected, minLength, maxLength, bounds, ProteinKmerGeometry.exact(5), positional, partitions, overlap, blacklist);
	}

	public static NegativePanel loadPackagedPanel(Path file, int expected, int minLength, int maxLength,
			FamilyTailCascade.ScalarBounds bounds, ProteinKmerGeometry geometry, boolean positional, int partitions, int overlap, DDLPentamerBlacklist blacklist){
		return loadPackagedPanel(file, expected, minLength, maxLength, bounds, geometry,
			positional ? ProteinPositionMode.FIXED : ProteinPositionMode.PLAIN, partitions, overlap, blacklist);
	}

	/** The query sketch coordinate rule is part of the immutable calibration panel identity. */
	public static NegativePanel loadPackagedPanel(Path file, int expected, int minLength, int maxLength,
			FamilyTailCascade.ScalarBounds bounds, ProteinKmerGeometry geometry, ProteinPositionMode positionMode,
			int partitions, int overlap, DDLPentamerBlacklist blacklist){
		return loadPackagedPanel(file, expected, minLength, maxLength, bounds, geometry, positionMode, partitions, overlap, blacklist, 1024);
	}

	/** Explicit width binds query sketches to the family rows before ranking or alignment. */
	public static NegativePanel loadPackagedPanel(Path file, int expected, int minLength, int maxLength,
			FamilyTailCascade.ScalarBounds bounds, ProteinKmerGeometry geometry, ProteinPositionMode positionMode,
			int partitions, int overlap, DDLPentamerBlacklist blacklist, int width){
		return loadPanel(file, expected, minLength, maxLength, bounds, geometry, positionMode, partitions, overlap, blacklist, false, width);
	}

	/** Also retains raw dimer counts for the separate legacy cosine gate, without recounting queries.
	 * Counts alone require 3200 bytes per query; callers must budget this optional storage and release
	 * the panel when both raw-cosine and packaged calibration are finished. Ordinary scoring omits it.
	 */
	public static NegativePanel loadPreparationPanel(Path file, int expected, int minLength, int maxLength,
			FamilyTailCascade.ScalarBounds bounds, boolean positional, int partitions, int overlap, DDLPentamerBlacklist blacklist){
		return loadPreparationPanel(file, expected, minLength, maxLength, bounds, ProteinKmerGeometry.exact(5), positional, partitions, overlap, blacklist);
	}

	public static NegativePanel loadPreparationPanel(Path file, int expected, int minLength, int maxLength,
			FamilyTailCascade.ScalarBounds bounds, ProteinKmerGeometry geometry, boolean positional, int partitions, int overlap, DDLPentamerBlacklist blacklist){
		return loadPreparationPanel(file, expected, minLength, maxLength, bounds, geometry,
			positional ? ProteinPositionMode.FIXED : ProteinPositionMode.PLAIN, partitions, overlap, blacklist);
	}

	/** Retains raw dimer counts as well as the explicitly identified query layout. */
	public static NegativePanel loadPreparationPanel(Path file, int expected, int minLength, int maxLength,
			FamilyTailCascade.ScalarBounds bounds, ProteinKmerGeometry geometry, ProteinPositionMode positionMode,
			int partitions, int overlap, DDLPentamerBlacklist blacklist){
		return loadPreparationPanel(file, expected, minLength, maxLength, bounds, geometry, positionMode, partitions, overlap, blacklist, 1024);
	}

	/** Retains raw counts without changing their arithmetic when the sketch width changes. */
	public static NegativePanel loadPreparationPanel(Path file, int expected, int minLength, int maxLength,
			FamilyTailCascade.ScalarBounds bounds, ProteinKmerGeometry geometry, ProteinPositionMode positionMode,
			int partitions, int overlap, DDLPentamerBlacklist blacklist, int width){
		return loadPanel(file, expected, minLength, maxLength, bounds, geometry, positionMode, partitions, overlap, blacklist, true, width);
	}

	private static NegativePanel loadPanel(Path file, int expected, int minLength, int maxLength,
			FamilyTailCascade.ScalarBounds bounds, ProteinKmerGeometry geometry, ProteinPositionMode positionMode, int partitions, int overlap,
			DDLPentamerBlacklist blacklist, boolean retainRaw, int width){
		if(width!=512 && width!=1024){throw new IllegalArgumentException("Negative panels require validated width512 or1024: "+width);}
		if(file==null || bounds==null || expected<1 || expected>MAX_PANEL_QUERIES){
			throw new IllegalArgumentException("Packaged negative panels require a file, scalar bounds and1.."+MAX_PANEL_QUERIES+" queries");
		}
		if(!Shared.AMINO_IN){throw new IllegalArgumentException("Configure amino-acid input before loading a packaged negative panel");}
		if(geometry==null || geometry.k()<4 || geometry.k()>6){throw new IllegalArgumentException("Packaged negative panels require base k from four through six");}
		if(positionMode==null){throw new IllegalArgumentException("Negative panels require an explicit position mode");}
		if(positionMode.positional()){DDLFamilyProfiles.validatePositionalLayout(partitions, overlap, geometry);}
		return new NegativePanel(loadNegatives(file, expected, minLength, maxLength, geometry, 0, width, bounds, positionMode, partitions, overlap, blacklist, retainRaw),
			bounds, geometry.name(), positionMode, partitions, overlap, blacklist==null ? "none" : blacklist.sha80(), retainRaw, width);
	}

	/** No query arrays escape; a panel can be shared by independent family workers. */
	public static final class NegativePanel{
		private NegativePanel(ArrayList<Negative> source, FamilyTailCascade.ScalarBounds bounds_, String geometry_,
				ProteinPositionMode positionMode_, int partitions_, int overlap_, String blacklistSha80_, boolean rawCounts_, int width_){
			assert(!source.isEmpty() && bounds_!=null) : "Combined25 family and query features require the same nonempty panel and bounds";
			queries=source.toArray(new Negative[source.size()]); bounds=bounds_;
			geometry=geometry_; positionMode=positionMode_; positional=positionMode.positional();
			partitions=partitions_; overlap=overlap_; blacklistSha80=blacklistSha80_;
			rawCounts=rawCounts_; width=width_;
			assert(width==512 || width==1024) : "Panel width must identify a validated query kernel before family scoring";
			for(Negative query:queries){
				assert(query.sketch.length==width) : "Index matching requires every query to use the panel's declared register width";
			}
			// Validate once before sharing: the exact-integer dot bound relies on nonnegative counts and exact totals.
			if(rawCounts){
				rawPairs=new long[queries.length][];
				for(int i=0; i<queries.length; i++){rawPairs[i]=rawPackedPairs(queries[i].spectrum, queries[i].norm);}
			}else{rawPairs=null;}
		}
		public int size(){return queries.length;}
		public int width(){return width;}
		public String queryId(int ordinal){return queries[ordinal].sequence.id;}
		/** Fills a caller-owned raw-DDL tile, widening unsigned short sketches only into scratch rows.
		 * Returned Query objects borrow the supplied scratch rows until all tile consumers finish.
		 */
		public void rawPairQueries(final int panel, final int from,
				final FamilyRawPairCollector.Query[] rows, final char[][] sketches){
			if(panel<0 || from<0 || rows==null || sketches==null || rows.length<1 || sketches.length<rows.length ||
					(long)from+rows.length>queries.length){
				throw new IllegalArgumentException("Invalid raw pair tile panel="+panel+" from="+from);
			}
			final HashSet<char[]> seen=new HashSet<char[]>(rows.length*2+1);
			for(int i=0; i<rows.length; i++){
				final char[] sketch=sketches[i];
				if(sketch==null || sketch.length!=width || !seen.add(sketch)){
					throw new IllegalArgumentException("Raw pair scratch rows must be distinct with width "+width);
				}
			}
			for(int i=0; i<rows.length; i++){
				final char[] sketch=sketches[i];
				final short[] source=queries[from+i].sketch;
				for(int j=0; j<width; j++){sketch[j]=(char)(source[j]&0xffff);}
				rows[i]=new FamilyRawPairCollector.Query(panel, from+i, sketch);
			}
		}
		/** Scans every query using original counts/norms, independently of sketch ranking or masks.
		 * No alignment or acceptance gate is applied. Concurrent calls may share this immutable panel;
		 * the caller must keep the prepared family's spectrum unchanged for the duration of each call.
		 */
		public double rawCosineMaximum(FamilyMemberPreparation.Prepared family){
			if(!rawCounts || family==null || family.members<1 || family.spectrum==null || family.spectrum.k!=2){
				throw new IllegalArgumentException("Raw cosine calibration requires a preparation panel and a live family's dimer counts");
			}
			family.spectrum.validate();
			final double norm=ProteinSpectrum.norm(family.spectrum);
			double maximum=0;
			for(int i=0; i<queries.length; i++){
				final Negative query=queries[i];
				final double cosine=rawCosineSparseValidated(query.spectrum, family.spectrum, query.norm, norm, rawPairs[i]);
				if(cosine>maximum){maximum=cosine;}
			}
			return maximum;
		}
		private final Negative[] queries;
		private final long[][] rawPairs;
		private final FamilyTailCascade.ScalarBounds bounds;
		private final String geometry;
		private final boolean positional, rawCounts;
		private final ProteinPositionMode positionMode;
		private final int partitions, overlap, width;
		private final String blacklistSha80;
	}

	/** Raw cosine for spectra whose validate() has succeeded; arrays and totals must remain immutable
	 * during scoring. Norms must come from those same arrays. The panel validates its private queries
	 * once, and rawCosineMaximum validates each borrowed family before scanning the panel.
	 * A sufficient total-product bound proves every integer product and partial sum exactly representable
	 * as a double; other inputs retain the original ordered floating-point calculation.
	 */
	public static double rawCosineValidated(ProteinSpectrum.Spectrum query, ProteinSpectrum.Spectrum family,
			double queryNorm, double familyNorm){
		if(query==null || family==null || query.k!=family.k || query.counts.length!=family.counts.length ||
				query.total<0 || family.total<0 || !Double.isFinite(queryNorm) || !Double.isFinite(familyNorm) ||
				queryNorm<0 || familyNorm<0){
			throw new IllegalArgumentException("Raw cosine requires matching validated spectra and finite nonnegative cached norms");
		}
		if(queryNorm==0 || familyNorm==0){return 0;}
		// Nonnegative terms imply dot<=sum(query)*sum(family); divide first to avoid bound overflow.
		if(family.total>0 && query.total>EXACT_DOUBLE_INTEGER_LIMIT/family.total){
			return ProteinSpectrum.cosine(query, family, queryNorm, familyNorm);
		}
		long dot=0;
		for(int i=0; i<query.counts.length; i++){dot+=query.counts[i]*family.counts[i];}
		return Math.min(1, (double)dot/queryNorm/familyNorm);
	}
	private static final long EXACT_DOUBLE_INTEGER_LIMIT=1L<<53;

	/** Fast raw-cosine path for immutable nonnegative 400-bin dimer spectra already validated by
	 * NegativePanel and rawCosineMaximum. Keeps the original zero-norm behavior, the original
	 * ordered floating fallback above the exact-integer total-product bound, and the original
	 * final division/max order while skipping query bins known to be zero.
	 */
	private static double rawCosineSparseValidated(ProteinSpectrum.Spectrum query, ProteinSpectrum.Spectrum family,
			double queryNorm, double familyNorm, long[] pairs){
		if(queryNorm==0 || familyNorm==0){return 0;}
		if(pairs==null || family.total>0 && query.total>EXACT_DOUBLE_INTEGER_LIMIT/family.total){
			return rawCosineValidated(query, family, queryNorm, familyNorm);
		}
		long dot=0;
		for(int i=0; i<pairs.length; i++){
			final long pair=pairs[i];
			dot+=(pair>>>9)*family.counts[(int)(pair&0x1FF)];
		}
		return Math.min(1, (double)dot/queryNorm/familyNorm);
	}

	/** Nonzero raw count/index pairs; a null row marks queries that need the dense fallback oracle. */
	private static long[] rawPackedPairs(ProteinSpectrum.Spectrum spectrum, double norm){
		if(spectrum==null || spectrum.k!=2 || !Double.isFinite(norm) || norm<0){
			throw new IllegalArgumentException("Raw query cache requires a validated dimer spectrum and finite norm");
		}
		spectrum.validate();
		if(spectrum.total>EXACT_DOUBLE_INTEGER_LIMIT){return null;}
		int nonzero=0;
		for(long count:spectrum.counts){if(count!=0){nonzero++;}}
		final long[] pairs=new long[nonzero];
		for(int i=0, j=0; i<spectrum.counts.length; i++){
			final long count=spectrum.counts[i];
			if(count!=0){pairs[j++]=(count<<9)|i;}
		}
		return pairs;
	}

	/** Shared fixed family features. The roster and sketch rows are copied; family payloads are
	 * borrowed and must remain immutable while workers run. Replacement rows must match the panel
	 * geometry/blacklist, as checked by the input loader or the enclosing preparation driver.
	 */
	public static final class PackagedCalibration{
		public PackagedCalibration(FamilyReference[] families_, char[][] replacement, NegativePanel panel_){
			this(families_, replacement, null, panel_);
		}
		/** Uses member-derived features directly; no provisional acceptance gates are invented.
		 * Prepared payloads are borrowed immutable. Their layout/blacklist must match the panel.
		 */
		public PackagedCalibration(FamilyMemberPreparation.Prepared[] prepared, NegativePanel panel){
			this(null, null, prepared, panel);
		}
		private PackagedCalibration(FamilyReference[] families, char[][] replacement,
				FamilyMemberPreparation.Prepared[] prepared, NegativePanel panel_){
			this(families, replacement, prepared, panel_, FamilyAlignmentMetric.LEGACY, null, null, null);
		}
		/** New-metric preparation uses the complete roster's pinned metadata with an explicit batch map.
		 * Hit floors are the deployed admission rule, not fitted here. Legacy entry points stay separate.
		 */
		public PackagedCalibration(FamilyMemberPreparation.Prepared[] prepared, NegativePanel panel,
				FamilyAlignmentMetric metric, FamilyHeuristicMetadata metadata, int[] ordinals, int[] hitFloors){
			this(null, null, prepared, panel, metric, metadata, ordinals, hitFloors);
		}
		private PackagedCalibration(FamilyReference[] families, char[][] replacement,
				FamilyMemberPreparation.Prepared[] prepared, NegativePanel panel_, FamilyAlignmentMetric metric_,
				FamilyHeuristicMetadata metadata, int[] ordinals, int[] hitFloors){
			if((families==null)==(prepared==null) || panel_==null){
				throw new IllegalArgumentException("Supply one live family representation and its negative panel");
			}
			final int count=prepared==null ? families.length : prepared.length;
			if(count==0 || (replacement!=null && replacement.length!=count)){
				throw new IllegalArgumentException("Packaged calibration requires live families, a panel and matching replacement rows");
			}
			if(metric_==null){throw new IllegalArgumentException("Negative calibration requires an explicit alignment contract");}
			metric=metric_;
			final boolean measured=!metric.same(FamilyAlignmentMetric.LEGACY);
			if(measured){
				if(prepared==null || metadata==null || !metadata.hasHydrophobicity() || ordinals==null || hitFloors==null ||
						ordinals.length!=count || hitFloors.length!=count){
					throw new IllegalArgumentException("Pinned-metric calibration requires prepared members, bound medians/KD and explicit batch admission");
				}
				// Metadata is already parsed and immutable; the enclosing run verifies its pins
				// before/after all batches rather than rereading whole tables for every batch.
				medians=new double[count]; hydrophobicities=new double[count]; minimumHits=hitFloors.clone();
			}else{
				if(metadata!=null || ordinals!=null || hitFloors!=null){throw new IllegalArgumentException("Legacy calibration must not silently ignore ranking inputs");}
				medians=hydrophobicities=null; minimumHits=null;
			}
			panel=panel_; consensuses=new ProteinSequence[count]; familyRanks=new int[count];
			final String[] names=new String[count]; final char[][] profiles=new char[count][];
			vectors=new double[count][25]; norms=new double[count]; final float[] mono=new float[20];
			for(int f=0; f<count; f++){
				final byte[] raw; final float[] spectrum;
				if(prepared==null){
					final FamilyReference family=families[f];
					if(family==null || family.members<1){throw new IllegalArgumentException("Empty families have no fitted negative upper: index="+f);}
					if(family.explicitAlignmentCalibration){throw new IllegalArgumentException("Legacy negative calibration cannot consume explicitly bound alignment states: "+family.id);}
					familyRanks[f]=family.rank;
					names[f]=family.id; profiles[f]=replacement==null ? family.ddl : replacement[f];
					consensuses[f]=family.consensus; raw=family.consensusRaw; spectrum=family.spectrum;
				}else{
					final FamilyMemberPreparation.Prepared family=prepared[f];
					if(family==null || family.members<1){throw new IllegalArgumentException("Empty prepared families have no fitted negative upper: index="+f);}
					familyRanks[f]=family.rank;
					if(!family.alignmentMetric.same(metric)){
						throw new IllegalArgumentException("Negative and member calibration require identical alignment metrics: "+family.id);
					}
					if(measured){
						final int ordinal=ordinals[f];
						if(minimumHits[f]<1 || minimumHits[f]>panel.width || metadata.rank(ordinal)!=family.rank ||
								!metadata.id(ordinal).equals(family.id) || metadata.members(ordinal)!=family.members){
							throw new IllegalArgumentException("Batch metadata and positive raw-hit floor must bind the exact prepared family: "+family.id);
						}
						medians[f]=metadata.median(ordinal); hydrophobicities[f]=metadata.hydrophobicity(ordinal);
					}
					if(family.width!=panel.width || !family.geometry.equals(panel.geometry) || !family.blacklistSha80.equals(panel.blacklistSha80) || (panel.positional &&
							(family.positionMode!=panel.positionMode || family.partitions!=panel.partitions || family.overlap!=panel.overlap))){
						throw new IllegalArgumentException("Prepared family and negative panel use different sketch width, layout or blacklist: "+family.id);
					}
					names[f]=family.id; profiles[f]=panel.positional ? family.positional : family.plain;
					consensuses[f]=family.consensus; raw=family.consensusRaw;
					spectrum=FamilyReferenceExport.frequencies(family.spectrum);
				}
				if(profiles[f]==null || profiles[f].length!=panel.width){
					throw new IllegalArgumentException("Packaged family profiles must match negative-panel width="+panel.width+": "+names[f]);
				}
				panel.bounds.combine(projectedMonomers(spectrum, mono), raw, false, vectors[f]);
				norms[f]=FamilyTailCascade.combinedNorm(vectors[f]);
			}
			index=new DDLFamilySearch.Index(names, profiles, null, false);
		}
		public PackagedWorker newWorker(int top){return new PackagedWorker(this, top);}
		public MetricWorker newMetricWorker(int top){return new MetricWorker(this, top);}
		private final ProteinSequence[] consensuses;
		private final NegativePanel panel;
		private final DDLFamilySearch.Index index;
		private final double[][] vectors;
		private final double[] norms;
		private final FamilyAlignmentMetric metric;
		private final double[] medians, hydrophobicities;
		private final int[] minimumHits;
		private final int[] familyRanks;
	}

	/** Synchronous selected-pair observation. Scores retain their original doubles; familyIndex is
	 * the supplied live-roster index, while queryOrdinal is the panel FASTA ordinal. No scratch escapes.
	 */
	public interface PairObserver{
		void observe(int familyIndex, int queryOrdinal, int ddl, double combinedCosine, double identityPct, double blosum);
	}

	/** One thread owns each worker. Every call scans the entire panel, selects topN with the
	 * existing ordinal tie rule, and aligns all selected pairs without acceptance/early-exit gates.
	 */
	public static final class PackagedWorker{
		private PackagedWorker(PackagedCalibration population_, int top_){
			if(!population_.metric.same(FamilyAlignmentMetric.LEGACY)){
				throw new IllegalArgumentException("Explicit alignment states require the metric-aware negative worker");
			}
			if(top_<1 || top_>population_.panel.size()){throw new IllegalArgumentException("Require1<=top<=panel size for exact family-side selection");}
			population=population_; top=top_; heap=new FamilyCascadeRanking.QueryScratch(top); hard=new int[top];
		}
		/** Returns detached maxima. Observer failure aborts this call; a later call resets the heap. */
		public FamilyResult scoreFamily(int f, PairObserver observer){
			if(f<0 || f>=population.consensuses.length){throw new IllegalArgumentException("Family index outside the supplied live roster: "+f);}
			final ProteinSequence reference=population.consensuses[f]; final Negative[] queries=population.panel.queries;
			final double[] vector=population.vectors[f]; final double norm=population.norms[f];
			heap.clear(); final long selectionStart=System.nanoTime(); double cosineMax=0;
			for(int q=0; q<queries.length; q++){
				final Negative query=queries[q]; final int ddl=population.index.scoreFamily(query.sketch, f);
				final double cosine=FamilyTailCascade.combinedCosine(query.combined, vector, query.combinedNorm, norm);
				cosineMax=Math.max(cosineMax, cosine); heap.offer(q, packagedRank(ddl, cosine), top);
			}
			final long selectedAt=System.nanoTime(); final int count=heap.drain(hard);
			if(count!=top){throw new IllegalStateException("Every live family must select exactly top negatives");}
			double identity=0, blosum=0;
			for(int i=0; i<count; i++){
				final int ordinal=hard[i]; final Negative q=queries[ordinal];
				FamilyCalibration.align(q.sequence, reference, q.selfScore, metrics);
				if(i==0 || metrics.identityPct>identity){identity=metrics.identityPct;} if(i==0 || metrics.blosum>blosum){blosum=metrics.blosum;}
				if(observer!=null){
					final int ddl=population.index.scoreFamily(q.sketch, f);
					final double cosine=FamilyTailCascade.combinedCosine(q.combined, vector, q.combinedNorm, norm);
					observer.observe(f, ordinal, ddl, cosine, metrics.identityPct, metrics.blosum);
				}
			}
			return new FamilyResult(f, queries.length, cosineMax, count, identity, blosum, selectedAt-selectionStart, 0, System.nanoTime()-selectedAt);
		}
		private final PackagedCalibration population;
		private final int top;
		private final FamilyCascadeRanking.QueryScratch heap;
		private final int[] hard;
		private final FamilyCalibration.Metrics metrics=new FamilyCalibration.Metrics();
	}

	/** Synchronous callback borrows the observation only until return; copy any retained fields. */
	public interface MetricPairObserver{
		void observe(int familyIndex, int queryOrdinal, int hits, double cosine, FamilyAlignmentMetric.Observation alignment);
	}

	/** One worker per thread, with reusable selection/alignment scratch and no per-pair allocation. */
	public static final class MetricWorker{
		private MetricWorker(PackagedCalibration population_, int top_){
			if(population_.metric.same(FamilyAlignmentMetric.LEGACY) || top_<1 || top_>population_.panel.size()){
				throw new IllegalArgumentException("Metric-aware calibration requires a new metric and1<=top<=panel size");
			}
			population=population_; top=top_; heap=new FamilyCascadeRanking.QueryScratch(top); hard=new int[top];
			alignment=population.metric.newWorker();
		}
		public MetricResult scoreFamily(int family, MetricPairObserver observer){
			return scoreFamily(family,observer,true);
		}
		/** Preparation exports alignment evidence, not the combined25 all-query maximum.
		 * Raw400 cosine calibration is a separate all-query scan and is not skipped here.
		 */
		public MetricResult scoreFamilyForPreparation(int family, MetricPairObserver observer){
			return scoreFamily(family,observer,false);
		}
		/** Aligns a sorted unique panel-local frontier without repeating all-query DDL selection.
		 * The returned scanned is the logical panel size for denominator parity; admitted and
		 * selected count this frontier, not every raw-positive candidate seen by the collector.
		 */
		public MetricResult scoreSelectedForPreparation(final int family, final int[] ordinals,
				final int count, final MetricPairObserver observer){
			checkFamily(family);
			final Negative[] queries=population.panel.queries;
			if(ordinals==null || count<0 || count>ordinals.length || count>queries.length){
				throw new IllegalArgumentException("Explicit negative frontier has invalid bounds");
			}
			Arrays.fill(counts, 0);
			double identity=0, blosum=0; int eligible=0;
			final long selectedAt=System.nanoTime();
			final double[] vector=observer==null ? null : population.vectors[family];
			final double norm=observer==null ? 0 : population.norms[family];
			int previous=-1;
			for(int i=0; i<count; i++){
				final int ordinal=ordinals[i];
				if(ordinal<=previous || ordinal>=queries.length){
					throw new IllegalArgumentException("Explicit negative ordinal out of panel bounds: "+ordinal);
				}
				previous=ordinal;
				final Negative query=queries[ordinal];
				alignment.measure(query.sequence, population.consensuses[family], query.selfScore, observation);
				counts[observation.state.ordinal()]++;
				if(observation.measured()){
					if(eligible==0 || observation.identityPct()>identity){identity=observation.identityPct();}
					if(eligible==0 || observation.blosum()>blosum){blosum=observation.blosum();}
					eligible++;
				}
				if(observer!=null){
					observer.observe(family, ordinal, population.index.scoreFamily(query.sketch, family),
						FamilyTailCascade.combinedCosine(query.combined, vector, query.combinedNorm, norm), observation);
				}
			}
			return new MetricResult(family, population.familyRanks[family], population.index.names[family], population.metric,
				queries.length, count, count, counts, 0, false, identity, blosum, 0, System.nanoTime()-selectedAt);
		}
		private MetricResult scoreFamily(int family, MetricPairObserver observer, boolean measureCosineMaximum){
			checkFamily(family);
			final Negative[] queries=population.panel.queries;
			final double[] vector=population.vectors[family]; final double norm=population.norms[family];
			heap.clear(); Arrays.fill(counts, 0);
			int admitted=0; double cosineMax=0;
			final long started=System.nanoTime();
			for(int q=0; q<queries.length; q++){
				final Negative query=queries[q]; final int hits=population.index.scoreFamily(query.sketch, family);
				if(!measureCosineMaximum && hits<population.minimumHits[family]){continue;}
				final double cosine=FamilyTailCascade.combinedCosine(query.combined, vector, query.combinedNorm, norm);
				if(measureCosineMaximum){cosineMax=Math.max(cosineMax, cosine);}
				if(hits<population.minimumHits[family]){continue;}
				if(!query.hydrophobicityDefined){throw new IllegalArgumentException("Admitted negative lacks measured query KD: "+query.sequence.id);}
				admitted++;
				final double similarity=FamilyHeuristicRanking.hydrophobicitySimilarity(query.hydrophobicity, population.hydrophobicities[family]);
				// Arm16 is the released fourth-root length +1.5 cosine +0.2 member-KD pivot.
				heap.offer(q, FamilyHeuristicRanking.score(16, hits, cosine, query.sequence.length(), population.medians[family], similarity), top);
			}
			final int selected=heap.drain(hard); final long selectedAt=System.nanoTime();
			if(selected!=Math.min(top, admitted)){throw new IllegalStateException("Every admitted negative must participate in deterministic top selection");}
			double identity=0, blosum=0; int eligible=0;
			for(int i=0; i<selected; i++){
				final int ordinal=hard[i]; final Negative query=queries[ordinal];
				alignment.measure(query.sequence, population.consensuses[family], query.selfScore, observation);
				counts[observation.state.ordinal()]++;
				if(observation.measured()){
					if(eligible==0 || observation.identityPct()>identity){identity=observation.identityPct();}
					if(eligible==0 || observation.blosum()>blosum){blosum=observation.blosum();}
					eligible++;
				}
				if(observer!=null){observer.observe(family, ordinal, population.index.scoreFamily(query.sketch, family),
					FamilyTailCascade.combinedCosine(query.combined, vector, query.combinedNorm, norm), observation);}
			}
			return new MetricResult(family, population.familyRanks[family], population.index.names[family], population.metric,
				queries.length, admitted, selected, counts, cosineMax,measureCosineMaximum,
				identity, blosum, selectedAt-started, System.nanoTime()-selectedAt);
		}
		/** Validates the family ordinal before using it in roster-indexed arrays. */
		private void checkFamily(final int family){
			if(family<0 || family>=population.consensuses.length){
				throw new IllegalArgumentException("Invalid negative-calibration family ordinal: "+family);
			}
		}
		private final PackagedCalibration population;
		private final int top;
		private final FamilyCascadeRanking.QueryScratch heap;
		private final int[] hard, counts=new int[FamilyAlignmentMetric.State.values().length];
		private final FamilyAlignmentMetric.Worker alignment;
		private final FamilyAlignmentMetric.Observation observation=new FamilyAlignmentMetric.Observation();
	}

	/** Detached negative evidence. Absence is explicit; normalized upper access requires measurements. */
	public static final class MetricResult{
		private MetricResult(int family_, int familyRank_, String familyId_, FamilyAlignmentMetric metric_, int scanned_, int admitted_, int selected_,
				int[] counts_, double cosineMax_, boolean cosineMaximumMeasured_, double identity_, double blosum_, long selectionNanos_, long alignmentNanos_){
			family=family_; metric=metric_; scanned=scanned_; admitted=admitted_; selected=selected_;
			familyRank=familyRank_; familyId=familyId_;
			counts=counts_.clone(); eligible=counts[FamilyAlignmentMetric.State.MEASURED.ordinal()];
			int total=0; for(int count:counts){if(count<0){throw new IllegalArgumentException("Negative alignment counts cannot be negative");} total=Math.addExact(total, count);}
			if(scanned<1 || admitted<0 || admitted>scanned || selected<0 || selected>admitted || total!=selected ||
					!Double.isFinite(cosineMax_) || cosineMax_<0 || cosineMax_>1 || !cosineMaximumMeasured_ && cosineMax_!=0 || !Double.isFinite(identity_) ||
					identity_<0 || identity_>100 || !Double.isFinite(blosum_) || eligible==0 && (identity_!=0 || blosum_!=0)){
				throw new IllegalArgumentException("Negative evidence must conserve selected outcomes and expose only finite measured maxima");
			}
			cosineMax=cosineMax_; identity=identity_; blosum=blosum_;
			cosineMaximumMeasured=cosineMaximumMeasured_;
			selectionNanos=selectionNanos_; alignmentNanos=alignmentNanos_;
		}
		public boolean measured(){return eligible>0;}
		public int ineligible(){return selected-eligible;}
		public int count(FamilyAlignmentMetric.State state){if(state==null){throw new IllegalArgumentException("Alignment count requires a state");} return counts[state.ordinal()];}
		public double identityUpper(){requireMeasured(); return identity;}
		public double blosumUpper(){requireMeasured(); return blosum;}
		public boolean hasCosineMaximum(){return cosineMaximumMeasured;}
		public double cosineMaximum(){
			if(!cosineMaximumMeasured){throw new IllegalStateException("Combined25 all-query cosine maximum was not measured for family "+family);}
			return cosineMax;
		}
		private void requireMeasured(){if(!measured()){throw new IllegalStateException("No eligible negative alignment upper exists for family "+family);}}
		public final int family, scanned, admitted, selected, eligible;
		public final int familyRank;
		public final String familyId;
		public final FamilyAlignmentMetric metric;
		/** Valid only when hasCosineMaximum(); preparation-only results store finite inactive zero.
		 * Existing scoreFamily callers always retain the full measured maximum.
		 */
		public final double cosineMax;
		private final boolean cosineMaximumMeasured;
		public final long selectionNanos, alignmentNanos;
		private final int[] counts;
		private final double identity, blosum;
	}

	private static float[] projectedMonomers(float[] spectrum,float[] scratch){
		if(spectrum.length==FamilyReferenceCodec.MONOMER_WIDTH){return spectrum;}
		if(spectrum.length==FamilyReferenceCodec.SPECTRUM_WIDTH){
			FamilyTailCascade.dimerRowSums(spectrum,scratch); return scratch;
		}
		throw new IllegalArgumentException("Family reference spectrum is neither projected nor dimer width");
	}
	public static double packagedRank(int ddl,double cosine){
		if(ddl<0 || ddl>1024 || !Double.isFinite(cosine) || cosine<0 || cosine>1){throw new IllegalArgumentException("Invalid family-side selection score");}
		return FamilyTailCascade.combinedRankingScore(ddl,cosine);
	}

	/** Records all DDL+cosine negative scores and retains top query ordinals per family. */
	public static final class NegativeSelector{
		public NegativeSelector(final int families, final int limit){
			if(families<1 || limit<1){throw new IllegalArgumentException("Invalid upper selector dimensions");}
			top=new FamilyCascadeRanking.QueryTopK(families, limit);
			cosineMax=new double[families]; observations=new int[families];
			Arrays.fill(cosineMax, 0);
		}
		public void clear(){
			top.clear();
			Arrays.fill(cosineMax, 0);
			Arrays.fill(observations, 0);
		}
		public void offer(final int family, final int ordinal, final int ddl, final double cosine){
			if(ddl<0 || !Double.isFinite(cosine) || cosine<0 || cosine>1){
				throw new IllegalArgumentException("Invalid negative score family="+family+" ordinal="+ordinal+
					" ddl="+ddl+" cosine="+cosine);
			}
			checkFamily(family);
			if(cosine>cosineMax[family]){cosineMax[family]=cosine;}
			observations[family]++;
			top.offer(family, ordinal, ddl+cosine);
		}
		public int drain(final int family, final int[] out){
			checkFamily(family);
			return top.drain(family, out);
		}
		public double cosineMax(final int family){
			checkFamily(family);
			if(observations[family]==0){throw new IllegalStateException("No negative observations for family: "+family);}
			return cosineMax[family];
		}
		public int observations(final int family){
			checkFamily(family);
			return observations[family];
		}
		public int families(){return observations.length;}
		public int limit(){return top.limit();}
		private void checkFamily(final int family){
			if(family<0 || family>=observations.length){throw new IllegalArgumentException("Invalid family index: "+family);}
		}
		private final FamilyCascadeRanking.QueryTopK top;
		private final double[] cosineMax;
		private final int[] observations;
	}

	/** Strict identity-percent upper gate;100 disables crossed gates without nonfinite state. */
	public static double enabledUpperGate(final double lowerGate, final double negativeMax){
		if(!Double.isFinite(lowerGate)){throw new IllegalArgumentException("Invalid lower gate: "+lowerGate);}
		if(!Double.isFinite(negativeMax)){throw new IllegalArgumentException("Invalid negative maximum: "+negativeMax);}
		return negativeMax>lowerGate ? negativeMax : FamilyGateCodec.DISABLED_UPPER;
	}

	private static FamilyResult[] scoreFamilies(DDLFamilySearch.Index index,Model[] models,
			ProteinSpectrum.Spectrum[] spectra,double[] norms,ArrayList<Negative> queries,
			int top,int requestedThreads,int shard,int shards){
		final FamilyResult[] results=new FamilyResult[index.families()];
		int selected=0;
		for(int family=0; family<results.length; family++){
			if(family%shards==shard){selected++;}
		}
		if(selected<1){throw new IllegalArgumentException("Selected family shard is empty: "+shard+"/"+shards);}
		final int threads=Math.min(requestedThreads, selected);
		final AtomicInteger next=new AtomicInteger(shard);
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			workers.add(new ProcessThread(next,index,models,spectra,norms,queries,results,top,shards,failure));
		}
		startAndJoin(workers,failure);
		if(failure.get()!=null){throw new RuntimeException("Negative upper worker failed", failure.get());}
		for(ProcessThread worker : workers){
			if(!worker.success){throw new IllegalStateException("Negative upper worker exited without success: "+worker.getName());}
		}
		for(int family=0; family<results.length; family++){
			final FamilyResult result=results[family];
			if(family%shards==shard){
				if(result==null || result.family!=family || result.observations!=queries.size()){
					throw new IllegalStateException("Selected family missing or incomplete after join: "+family);
				}
			}else if(result!=null){
				throw new IllegalStateException("Unselected family was processed: "+family);
			}
		}
		return results;
	}

	/** Start all workers and always join the successfully-started prefix. */
	static void startAndJoin(ArrayList<? extends Thread> workers,AtomicReference<Throwable> failure){
		final ArrayList<Thread> started=new ArrayList<Thread>(workers.size());
		try{
			for(Thread worker : workers){
				worker.start();
				started.add(worker);
			}
		}catch(Throwable t){failure.compareAndSet(null, t);}
		finally{ThreadWaiter.waitForThreadsToFinish(started);}
	}

	private static DDLFamilySearch.Index loadIndex(ReferenceFiles references,ProteinKmerGeometry geometry,long seed,int width,Path work) throws Exception{
		final Path list=Files.createTempFile(work,"family_negative_upper_profiles",".tsv");
		try{
			write(list,references.profiles.toString()+"\t"+references.profilesHash+"\t"+references.summaryHash+"\n");
			return DDLFamilySearch.loadForSearch(list,width,"nonempty",geometry,seed,false,false,false);
		}finally{Files.deleteIfExists(list);}
	}

	private static LinkedSpectra loadSpectra(Path file,String[] names) throws Exception{
		final HashMap<String,ProteinSpectrum.Spectrum> loaded=ProteinSpectrum.load(file.toString());
		final ProteinSpectrum.Spectrum[] spectra=new ProteinSpectrum.Spectrum[names.length];
		for(int i=0; i<names.length; i++){
			spectra[i]=loaded.remove(names[i]);
			if(spectra[i]==null || spectra[i].k!=2){
				throw new IllegalArgumentException("Missing or invalid dimer spectrum: "+names[i]);
			}
		}
		if(!loaded.isEmpty()){
			throw new IllegalArgumentException("Unmatched spectra outside DDL profile roster: "+loaded.keySet());
		}
		return new LinkedSpectra(spectra);
	}

	private static Models loadModels(Path file,DDLFamilySearch.Index index,int maxReference) throws Exception{
		final Model[] models=new Model[index.families()];
		final HashMap<String,Integer> storage=new HashMap<String,Integer>(index.families()*2);
		for(int i=0; i<index.names.length; i++){
			if(!index.names[i].startsWith("tail:") || storage.put(index.names[i].substring(5),i)!=null){
				throw new IllegalArgumentException("Duplicate or non-tail DDL family: "+index.names[i]);
			}
		}
		final ByteFile input=ByteFile.makeByteFile1(file.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		int rank=0,found=0;
		try{
			header(input,HbmTailLibrary.MODELS_HEADER);
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line);
				if(lp.terms()!=9 || lp.parseInt(0)!=rank){
					throw new IllegalArgumentException("Fixed model roster is not canonical");
				}
				final String id=lp.parseString(1);
				final Integer indexSlot=storage.remove(id);
				if(indexSlot!=null){
					final Path consensus=Paths.get(lp.parseString(2)).resolve("consensus.faa");
					final String hash=lp.parseString(6);
					pin(consensus,hash);
					models[indexSlot.intValue()]=new Model(rank,id,consensus(consensus,maxReference,id));
					pin(consensus,hash);
					found++;
				}
				rank++;
			}
		}finally{if(input.close()){throw new IllegalStateException("Fixed models read failed");}}
		if(!storage.isEmpty() || found!=models.length){
			throw new IllegalArgumentException("Fixed models missing profile families: "+storage.keySet());
		}
		return new Models(models);
	}

	private static Models loadModels(Path file,Path consensus,DDLFamilySearch.Index index,int maxReference) throws Exception{
		final Model[] models=new Model[index.families()];
		final int[] ranks=new int[index.families()];
		Arrays.fill(ranks,-1);
		final String[] ids=new String[index.families()];
		final HashMap<String,Integer> storage=new HashMap<String,Integer>(index.families()*2);
		final HashSet<String> allIds=new HashSet<String>();
		for(int i=0; i<index.names.length; i++){
			if(!index.names[i].startsWith("tail:") || storage.put(index.names[i].substring(5),i)!=null){
				throw new IllegalArgumentException("Duplicate or non-tail DDL family: "+index.names[i]);
			}
		}
		final ByteFile input=ByteFile.makeByteFile1(file.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		int rank=0,found=0;
		try{
			header(input,HbmTailLibrary.MODELS_HEADER);
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line);
				if(lp.terms()!=9 || lp.parseInt(0)!=rank){throw new IllegalArgumentException("Fixed model roster is not canonical");}
				final String id=lp.parseString(1);
				if(!allIds.add(id)){throw new IllegalArgumentException("Duplicate fixed model family: "+id);}
				final Integer indexSlot=storage.remove(id);
				if(indexSlot!=null){
					ranks[indexSlot.intValue()]=rank; ids[indexSlot.intValue()]=id; found++;
				}
				rank++;
			}
		}finally{if(input.close()){throw new IllegalStateException("Fixed models read failed");}}
		if(!storage.isEmpty() || found!=models.length){
			throw new IllegalArgumentException("Fixed models missing profile families: "+storage.keySet());
		}
		final HashMap<String,ProteinSequence> sequences=loadCombinedConsensus(consensus,maxReference,allIds);
		for(int i=0; i<ids.length; i++){
			final ProteinSequence sequence=sequences.remove(ids[i]);
			if(sequence==null){throw new IllegalArgumentException("Missing combined consensus: "+ids[i]);}
			models[i]=new Model(ranks[i],ids[i],sequence);
		}
		return new Models(models);
	}

	private static HashMap<String,ProteinSequence> loadCombinedConsensus(Path path,int cap,HashSet<String> expected){
		final HashMap<String,ProteinSequence> map=new HashMap<String,ProteinSequence>();
		final FileFormat format=FileFormat.testInput(path.toString(),FileFormat.FASTA,null,false,true);
		final Streamer input=StreamerFactory.makeStreamer(format,null,true,-1,false,true,1);
		input.start();
		try{
			for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
				for(Read read : batch.list){
					if(read==null || read.bases==null || read.mate!=null || read.bases.length<1 || read.bases.length>cap){
						throw new IllegalArgumentException("Invalid combined consensus record");
					}
					final String id=firstToken(read.id);
					if(!expected.remove(id)){throw new IllegalArgumentException("Duplicate or foreign combined consensus: "+id);}
					map.put(id,new ProteinSequence(id,read.bases));
				}
			}
		}finally{if(ReadWrite.closeStream(input)){throw new IllegalStateException("Combined consensus read failed");}}
		if(!expected.isEmpty()){throw new IllegalArgumentException("Missing combined consensus: "+expected);}
		return map;
	}

	private static ProteinSequence consensus(Path path,int cap,String family){
		final FileFormat format=FileFormat.testInput(path.toString(),FileFormat.FASTA,null,false,true);
		final Streamer input=StreamerFactory.makeStreamer(format,null,true,-1,false,true,1);
		input.start();
		ProteinSequence found=null;
		try{
			for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
				for(Read read : batch.list){
					if(read==null || read.bases==null || read.mate!=null || read.bases.length<1 ||
							read.bases.length>cap || found!=null){
						throw new IllegalArgumentException("Fixed consensus must contain exactly one bounded protein");
					}
					final String id=firstToken(read.id);
					if(!id.equals(family)){
						throw new IllegalArgumentException("Fixed consensus ID "+id+" differs from models family "+family);
					}
					found=new ProteinSequence(id,read.bases);
				}
			}
		}finally{if(ReadWrite.closeStream(input)){throw new IllegalStateException("Fixed consensus read failed");}}
		if(found==null){throw new IllegalArgumentException("Missing fixed consensus: "+path);}
		return found;
	}

	private static ArrayList<Negative> loadNegatives(Path file,int expected,int minLength,int maxLength,int k,long seed,int width){
		return loadNegatives(file,expected,minLength,maxLength,k,seed,width,null,false,1,0,null);
	}

	private static void validateQueryIds(Path file,NegativePanel queries){
		final ByteFile input=ByteFile.makeByteFile1(file.toString(),false);
		final HashSet<String> seen=new HashSet<String>(queries.size()*2);
		int row=0;
		try{
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				if(row>=queries.size()){throw new IllegalArgumentException("Too many query IDs");}
				final String id=queries.queryId(row);
				if(!equals(line,id)){throw new IllegalArgumentException("Query ID mismatch at row "+row+": "+id);}
				if(!seen.add(id)){throw new IllegalArgumentException("Duplicate query ID: "+id);}
				row++;
			}
		}finally{if(input.close()){throw new IllegalStateException("Query ID list read failed");}}
		if(row!=queries.size()){throw new IllegalArgumentException("Wrong query ID count: "+row);}
	}
	private static ArrayList<Negative> loadNegatives(Path file,int expected,int minLength,int maxLength,int k,long seed,int width,FamilyTailCascade.ScalarBounds bounds){
		return loadNegatives(file,expected,minLength,maxLength,k,seed,width,bounds,false,1,0,null);
	}
	private static ArrayList<Negative> loadNegatives(Path file,int expected,int minLength,int maxLength,int k,long seed,int width,
			FamilyTailCascade.ScalarBounds bounds,boolean positional,int partitions,int overlap,DDLPentamerBlacklist blacklist){
		return loadNegatives(file, expected, minLength, maxLength, ProteinKmerGeometry.exact(k), seed, width, bounds, positional, partitions, overlap, blacklist, false);
	}
	private static ArrayList<Negative> loadNegatives(Path file,int expected,int minLength,int maxLength,int k,long seed,int width,
			FamilyTailCascade.ScalarBounds bounds,boolean positional,int partitions,int overlap,DDLPentamerBlacklist blacklist,boolean retainRaw){
		return loadNegatives(file, expected, minLength, maxLength, ProteinKmerGeometry.exact(k), seed, width, bounds, positional, partitions, overlap, blacklist, retainRaw);
	}
	private static ArrayList<Negative> loadNegatives(Path file,int expected,int minLength,int maxLength,ProteinKmerGeometry geometry,long seed,int width,
			FamilyTailCascade.ScalarBounds bounds,boolean positional,int partitions,int overlap,DDLPentamerBlacklist blacklist,boolean retainRaw){
		return loadNegatives(file, expected, minLength, maxLength, geometry, seed, width, bounds,
			positional ? ProteinPositionMode.FIXED : ProteinPositionMode.PLAIN, partitions, overlap, blacklist, retainRaw);
	}
	private static ArrayList<Negative> loadNegatives(Path file,int expected,int minLength,int maxLength,ProteinKmerGeometry geometry,long seed,int width,
			FamilyTailCascade.ScalarBounds bounds,ProteinPositionMode positionMode,int partitions,int overlap,DDLPentamerBlacklist blacklist,boolean retainRaw){
		if(expected<1 || minLength<1 || maxLength<minLength){
			throw new IllegalArgumentException("Invalid negative control bounds");
		}
		final ArrayList<Negative> queries=new ArrayList<Negative>();
		final char[] sketch=new char[width];
		final double[] scalarScratch=bounds==null ? null : new double[ConsensusScalars.SCALAR_COUNT];
		final FileFormat format=FileFormat.testInput(file.toString(),FileFormat.FASTA,null,false,true);
		final Streamer input=StreamerFactory.makeStreamer(format,null,true,-1,false,true,1);
		input.start();
		try{
			for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
				for(Read read : batch.list){
					if(read==null || read.bases==null || read.mate!=null){
						throw new IllegalArgumentException("Negative controls must be unpaired proteins");
					}
					if(queries.size()>=expected){throw new IllegalArgumentException("Too many negative controls");}
					final ProteinSequence seq=new ProteinSequence(firstToken(read.id),read.bases);
					if(seq.length()<minLength || seq.length()>maxLength){
						throw new IllegalArgumentException("Out-of-bounds negative length "+seq.length()+": "+seq.id);
					}
					final byte[] raw=Arrays.copyOf(read.bases,seq.length());
					final int self=ProteinSearcher.selfScore(seq);
					if(self<=0){throw new IllegalArgumentException("Negative query has nonpositive self score: "+seq.id);}
					final ProteinSpectrum.Spectrum spectrum=new ProteinSpectrum.Spectrum(2);
					spectrum.add(raw);
					if(positionMode.positional()){
						DDLFamilyProfiles.sketchPositional(raw,geometry,seed,partitions,overlap,positionMode,sketch,blacklist);
					}else{DDLFamilyProfiles.sketch(raw,geometry,seed,sketch,blacklist);}
					final short[] stored=new short[sketch.length];
					for(int i=0; i<sketch.length; i++){stored[i]=(short)sketch[i];}
					if(bounds==null){queries.add(new Negative(raw,seq,self,stored,spectrum,ProteinSpectrum.norm(spectrum),null));}
					else{
						final float[] mono=new float[20]; FamilyTailCascade.dimerRowSums(FamilyReferenceExport.frequencies(spectrum),mono);
						final int scalarMask=ConsensusScalars.computeScalarMask(raw, raw.length, scalarScratch, 0);
						final double[] combined=new double[25]; bounds.combineComputed(mono, true, scalarScratch, scalarMask, combined);
						queries.add(new Negative(raw, seq, self, stored, retainRaw ? spectrum : null,
							retainRaw ? ProteinSpectrum.norm(spectrum) : 0, combined,
							(scalarMask&(1<<ConsensusScalars.HYDRO))!=0, scalarScratch[ConsensusScalars.HYDRO]));
					}
				}
			}
		}finally{if(ReadWrite.closeStream(input)){throw new IllegalStateException("Negative FASTA read failed");}}
		if(queries.size()!=expected){
			throw new IllegalArgumentException("Wrong negative query count: "+queries.size());
		}
		return queries;
	}

	private static ReferenceFiles loadReferenceManifest(Path manifest) throws Exception{
		final Path dir=manifest.getParent();
		final java.util.List<String> lines=Files.readAllLines(manifest,StandardCharsets.US_ASCII);
		if(lines.size()!=8 || !lines.get(0).equals("file\tsha80")){
			throw new IllegalArgumentException("Require exact full-reference result manifest");
		}
		final HashMap<String,Path> files=new HashMap<String,Path>();
		final HashMap<String,String> hashes=new HashMap<String,String>();
		for(int i=1;i<lines.size();i++){
			final String[] fields=lines.get(i).split("\t",-1);
			if(fields.length!=2){throw new IllegalArgumentException("Malformed reference manifest row");}
			final String logical=logicalReferenceName(fields[0]);
			if(files.put(logical,dir.resolve(fields[0]))!=null || hashes.put(logical,fields[1])!=null){
				throw new IllegalArgumentException("Duplicate reference manifest logical file: "+logical);
			}
			pin(dir.resolve(fields[0]),fields[1]);
		}
		for(String name : REFERENCE_FILES){
			if(!files.containsKey(name)){throw new IllegalArgumentException("Missing reference file: "+name);}
		}
		if(!new String(Files.readAllBytes(files.get("PASS")),StandardCharsets.US_ASCII).equals("FAMILY_FULL_REFERENCES_PASS\n")){
			throw new IllegalArgumentException("Full-member references have not passed");
		}
		return new ReferenceFiles(files.get("profiles.tsv"),hashes.get("profiles.tsv"),
			files.get("spectra.tsv"),hashes.get("summary.tsv"));
	}

	private static String logicalReferenceName(String name){
		if(name.equals("profiles.tsv.gz")){return "profiles.tsv";}
		if(name.equals("spectra.tsv.gz")){return "spectra.tsv";}
		for(String allowed : REFERENCE_FILES){if(name.equals(allowed)){return name;}}
		throw new IllegalArgumentException("Unexpected reference manifest file: "+name);
	}

	private static void header(ByteFile input,String value){
		if(!Arrays.equals(input.nextLine(),value.getBytes(StandardCharsets.US_ASCII))){
			throw new IllegalArgumentException("Unexpected input header: "+value);
		}
	}

	private static Path input(HashMap<String,String> options,String key) throws Exception{
		final Path p=Paths.get(required(options,key)).toAbsolutePath().normalize();
		pin(p,required(options,key+"sha80"));
		return p;
	}

	private static void pin(Path path,String hash) throws Exception{
		if(!Files.isRegularFile(path) || Files.isSymbolicLink(path)){
			throw new IllegalArgumentException("Require a regular input: "+path);
		}
		if(!DigestSuffix.file(path.toString()).equals(DigestSuffix.requireSuffix(hash,"input sha80"))){
			throw new IllegalArgumentException("Input pin mismatch: "+path);
		}
	}

	private static String firstToken(String id){
		if(id==null || id.isEmpty()){throw new IllegalArgumentException("Missing FASTA ID");}
		int end=0; while(end<id.length() && !Character.isWhitespace(id.charAt(end))){end++;}
		return id.substring(0,end);
	}

	private static boolean equals(byte[] line,String s){
		if(line.length!=s.length()){return false;}
		for(int i=0;i<line.length;i++){if(line[i]!=(byte)s.charAt(i)){return false;}}
		return true;
	}

	private static long number(HashMap<String,String> o,String key){
		return Long.parseLong(required(o,key));
	}

	private static int integer(HashMap<String,String> o,String key){
		return Math.toIntExact(number(o,key));
	}

	private static int integer(HashMap<String,String> o,String key,int defaultValue){
		final String value=o.get(key);
		return value==null ? defaultValue : Math.toIntExact(Long.parseLong(value));
	}

	private static int threads(HashMap<String,String> o){
		final String t=o.get("t"), threads=o.get("threads");
		if(t!=null && threads!=null && !t.equals(threads)){throw new IllegalArgumentException("t and threads must agree");}
		final String value=t==null ? threads : t;
		final int result=value==null ? 1 : Math.toIntExact(Long.parseLong(value));
		if(result<1 || result>64){throw new IllegalArgumentException("threads must be in 1..64: "+result);}
		return result;
	}

	private static boolean optionalBoolean(HashMap<String,String> o,String key,boolean defaultValue){
		final String value=o.get(key);
		return value==null ? defaultValue : Parse.parseBoolean(value);
	}

	private static boolean containsAny(HashMap<String,String> options,String[] keys){
		for(String key : keys){if(options.containsKey(key)){return true;}}
		return false;
	}

	private static String required(HashMap<String,String> o,String key){
		final String value=o.get(key);
		if(value==null || value.isEmpty()){throw new IllegalArgumentException("Missing option: "+key);}
		return value;
	}

	private static ByteStreamWriter writer(Path path){
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(),false,false,false);
		writer.start();
		return writer;
	}

	private static void write(Path path,String value) throws Exception{
		Files.write(path,value.getBytes(StandardCharsets.US_ASCII));
	}

	private static void writeTimings(Path out,FamilyResult[] results,Model[] models,int families,int negatives,
			long cosineCalls,long alignments,int top,int threads,int shard,int shards,long started,long loaded) throws Exception{
		final StringBuilder table=new StringBuilder("rank\tfamily_id\tddl_score_seconds\tcosine_seconds\tforced_alignment_seconds\n");
		long ddl=0,cosine=0,aligned=0;
		int processed=0;
		for(FamilyResult result : results){
			if(result==null){continue;}
			ddl+=result.ddlNanos; cosine+=result.cosineNanos; aligned+=result.alignNanos; processed++;
			table.append(models[result.family].rank).append('\t').append(models[result.family].id).append('\t')
				.append(result.ddlNanos*1e-9).append('\t')
				.append(result.cosineNanos*1e-9).append('\t').append(result.alignNanos*1e-9).append('\n');
		}
		write(out.resolve("timing.tsv"),table.toString());
		write(out.resolve("summary.tsv"),"families\tprocessed_families\tnegatives\tcosine_calls\talignments\ttop\tthreads\tshard\tshards\tload_seconds\tddl_score_seconds\tcosine_seconds\talign_seconds\ttotal_seconds\n"+
			families+"\t"+processed+"\t"+negatives+"\t"+cosineCalls+"\t"+alignments+"\t"+top+"\t"+
			threads+"\t"+shard+"\t"+shards+"\t"+((loaded-started)*1e-9)+"\t"+
			(ddl*1e-9)+"\t"+(cosine*1e-9)+"\t"+(aligned*1e-9)+"\t"+((System.nanoTime()-started)*1e-9)+"\n");
	}

	private static final class ReferenceFiles{
		ReferenceFiles(Path profiles_,String profilesHash_,Path spectra_,String summaryHash_){
			profiles=profiles_; profilesHash=profilesHash_; spectra=spectra_; summaryHash=summaryHash_;
		}
		final Path profiles,spectra;
		final String profilesHash,summaryHash;
	}

	private static final class LinkedSpectra{
		LinkedSpectra(ProteinSpectrum.Spectrum[] values_){values=values_;}
		final ProteinSpectrum.Spectrum[] values;
	}

	private static final class Models{
		Models(Model[] values_){values=values_;}
		final Model[] values;
	}

	private static final class ProcessThread extends Thread{
		ProcessThread(AtomicInteger next_,DDLFamilySearch.Index index_,Model[] models_,
				ProteinSpectrum.Spectrum[] spectra_,double[] norms_,ArrayList<Negative> queries_,
				FamilyResult[] results_,int top_,int shards_,AtomicReference<Throwable> failure_){
			next=next_; index=index_; models=models_; spectra=spectra_; norms=norms_;
			queries=queries_; results=results_; top=top_; shards=shards_; failure=failure_;
			setName("family-negative-upper-"+next_.get()+"-"+System.identityHashCode(this));
			setDaemon(true);
		}
		@Override
		public void run(){
			try{
				final FamilyCascadeRanking.QueryScratch scratch=new FamilyCascadeRanking.QueryScratch(top);
				final FamilyCalibration.Metrics metrics=new FamilyCalibration.Metrics();
				final int[] hard=new int[top];
				for(int family=next.getAndAdd(shards); family<results.length && failure.get()==null;
						family=next.getAndAdd(shards)){
					double cosineMax=0;
					long ddlNanos=0,cosineNanos=0,alignNanos=0;
					scratch.clear();
					for(int ordinal=0; ordinal<queries.size(); ordinal++){
						final Negative query=queries.get(ordinal);
						final long ddlStart=System.nanoTime();
						final int ddl=index.scoreFamily(query.sketch,family);
						final long cosineStart=System.nanoTime();
						final double cosine=ProteinSpectrum.cosine(query.spectrum,spectra[family],query.norm,norms[family]);
						final long done=System.nanoTime();
						ddlNanos+=cosineStart-ddlStart; cosineNanos+=done-cosineStart;
						if(cosine>cosineMax){cosineMax=cosine;}
						scratch.offer(ordinal, ddl+cosine, top);
					}
					double identity=0, blosum=0;
					final int aligned=scratch.drain(hard);
					if(aligned<1){throw new IllegalStateException("No selected negative alignment for family: "+family);}
					final long alignStart=System.nanoTime();
					for(int i=0; i<aligned; i++){
						final Negative query=queries.get(hard[i]);
						FamilyCalibration.align(query.sequence,models[family].sequence,query.selfScore,metrics);
						if(i==0 || metrics.identityPct>identity){identity=metrics.identityPct;}
						if(i==0 || metrics.blosum>blosum){blosum=metrics.blosum;}
					}
					alignNanos=System.nanoTime()-alignStart;
					results[family]=new FamilyResult(family,queries.size(),cosineMax,aligned,identity,blosum,
						ddlNanos,cosineNanos,alignNanos);
				}
				success=true;
			}catch(Throwable e){failure.compareAndSet(null, e);}
		}
		private final AtomicInteger next;
		private final DDLFamilySearch.Index index;
		private final Model[] models;
		private final ProteinSpectrum.Spectrum[] spectra;
		private final double[] norms;
		private final ArrayList<Negative> queries;
		private final FamilyResult[] results;
		private final int top,shards;
		private final AtomicReference<Throwable> failure;
		private boolean success;
	}

	/** Detached maxima/counts; family is the live-roster index, and timings are nanoseconds.
	 * Packaged results use combined25 cosine, legacy results use raw dimer cosine. These are not interchangeable.
	 */
	public static final class FamilyResult{
		private FamilyResult(int family_,int observations_,double cosineMax_,int aligned_,
				double identity_,double blosum_,long ddlNanos_,long cosineNanos_,long alignNanos_){
			if(observations_<1 || aligned_<1 || !Double.isFinite(cosineMax_) || cosineMax_<0 || cosineMax_>1 ||
					!Double.isFinite(identity_) || identity_<0 || identity_>100 || !Double.isFinite(blosum_)){
				throw new IllegalArgumentException("Upper maxima require finite observed scores: family="+family_+" observations="+observations_+" aligned="+aligned_);
			}
			family=family_; observations=observations_; cosineMax=cosineMax_; aligned=aligned_;
			identity=identity_; blosum=blosum_; ddlNanos=ddlNanos_; cosineNanos=cosineNanos_;
			alignNanos=alignNanos_;
		}
		public final int family,observations,aligned;
		public final double cosineMax,identity,blosum;
		public final long ddlNanos,cosineNanos,alignNanos;
	}

	private static final class Model{
		Model(int rank_,String id_,ProteinSequence sequence_){rank=rank_; id=id_; sequence=sequence_;}
		final int rank;
		final String id;
		final ProteinSequence sequence;
	}

	private static final class Negative{
		Negative(byte[] raw_,ProteinSequence sequence_,int selfScore_,short[] sketch_,
				ProteinSpectrum.Spectrum spectrum_,double norm_,double[] combined_){
			this(raw_, sequence_, selfScore_, sketch_, spectrum_, norm_, combined_, false, 0);
		}
		Negative(byte[] raw_,ProteinSequence sequence_,int selfScore_,short[] sketch_,
				ProteinSpectrum.Spectrum spectrum_,double norm_,double[] combined_,boolean hydrophobicityDefined_,double hydrophobicity_){
			raw=raw_; sequence=sequence_; selfScore=selfScore_;
			sketch=sketch_; spectrum=spectrum_; norm=norm_;
			combined=combined_; combinedNorm=combined==null ? 0 : FamilyTailCascade.combinedNorm(combined);
			hydrophobicityDefined=hydrophobicityDefined_; hydrophobicity=hydrophobicity_;
		}
		final byte[] raw;
		final ProteinSequence sequence;
		final int selfScore;
		final short[] sketch;
		final ProteinSpectrum.Spectrum spectrum;
		final double norm;
		final double[] combined;
		final double combinedNorm;
		final boolean hydrophobicityDefined;
		final double hydrophobicity;
	}

	/** Resource ceiling for a complete in-memory packaged panel, not a sampling limit.
	 * Query ordinals/counts and validateQueryIds' doubled capacity remain within int;
	 * comparison totals use long and worker heaps depend on top, not panel size.
	 * At this ceiling, short[1024] sketches alone use4.096GB; optional long[400]
	 * raw dimer counts add6.4GB. Callers must also budget sequences, vectors, IDs,
	 * object overhead and family state. This ceiling does not guarantee a heap fit.
	 */
	public static final int MAX_PANEL_QUERIES=2000000;

	private static final String[] REFERENCE_FILES={"profiles.tsv","spectra.tsv","members.tsv","summary.tsv","recipe.txt","timing.tsv","PASS"};
	private static final java.util.List<String> OPTIONS=Arrays.asList("references","referencessha80","models","modelssha80",
		"negatives","negativessha80","out","top","k","seed","width","maxreads","maxreferenceaa","minnegativeaa","maxnegativeaa",
		"t","threads","shard","shards","familyreferences","familyreferencessha80","maxfamilies","queryprefix","queryids","queryidssha80","scalarbounds","scalarboundssha80","scalarscale",
		"aabits","aapartition","ddlprofiles","ddlprofilessha80","positional","positionmode","partitions","overlap","blacklist","blacklistsha80","consensus","consensussha80");
	private static final String[] PACKAGED_REPLACEMENT_OPTIONS={"ddlprofiles","ddlprofilessha80","positional","positionmode","partitions","overlap","blacklist","blacklistsha80"};
	public static final String PACKAGED_RAW_HEADER="roster_rank\tfamily_id\tquery_ordinal\tquery_id\tddl_matches\tcombined_cosine\tidentity_pct\tblosum_a_sq";
}
