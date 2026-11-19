package prot;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import ddl.DDLFamilyHoldout;
import ddl.DDLFamilySearch.StructuralInputs;
import ddl.DDLFamilyProfiles;
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

/** All-member family lower cutoffs against fixed truth-label consensuses.
 * No holdout filtering is performed; full-member calibration is not independent validation.
 * One validated bucket load per mode; independent family workers, ordered output, no search.
 * @author Sayu, Keqing
 */
public final class FamilyCalibration{

	private FamilyCalibration(){}

	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args,FamilyCalibration.class,false);
			try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	public static void run(String[] args) throws Exception{
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			require(pair.length==2,"Require flag=value arguments");
			final String parsed=pair[0].toLowerCase(Locale.ROOT), key=parsed.equals("threads") ? "t" : parsed;
			require(OPTIONS.contains(key) && options.put(key,pair[1])==null,"Unknown or duplicate option: "+key);
		}
		if(!options.containsKey("t")){options.put("t","1");}
		final MemberLengthRange memberRange=MemberLengthRange.from(options);
		final int threads=integer(options,"t"); require(threads>0,"Family worker count must be positive");
		final int nf=integer(options,"expectedfamilies"), batches=integer(options,"expectedbatches");
		final int bucket=integer(options,"bucket"), buckets=integer(options,"buckets");
		final String mode=required(options,"mode");
		require(mode.equals("references") || mode.equals("lower"),"Require mode=references or lower");
		final boolean references=mode.equals("references");
		require(options.containsKey("blacklist")==options.containsKey("blacklistsha80"),"Blacklist path and fingerprint are required together");
		require(references || !options.containsKey("blacklist"),"Pentamer blacklist applies only to DDL reference construction");
		final boolean retainScores=options.containsKey("retainscores") && Parse.parseBoolean(options.get("retainscores"));
		require(!references || !retainScores,"retainscores applies only to member lower calibration");
		if(references){for(String key : new String[]{"models","modelssha80","references","referencessha80","maxreferenceaa"}){
			require(!options.containsKey(key),"Lower-only option supplied to reference construction: "+key);
		}}
		final int maxReference=references ? 200 : integer(options,"maxreferenceaa");
		final ProteinPositionMode positionMode=ProteinPositionMode.resolve(options.get("positionmode"), options.get("sketch"),
			options.containsKey("positional") ? Boolean.valueOf(Parse.parseBoolean(options.get("positional"))) : null, ProteinPositionMode.PLAIN);
		final String sketch=positionMode.sketch;
		require(!(options.containsKey("aabits") && options.containsKey("aapartition")),"Use either aabits or aapartition, not both");
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(options.containsKey("k") ? integer(options,"k") : 5,
			options.containsKey("aabits") ? options.get("aabits") : options.get("aapartition"));
		require(references || sketch.equals(DDLFamilyProfiles.SKETCH_PLAIN),"Lower scoring requires the plain reference recipe");
		final int partitions=options.containsKey("partitions") ? integer(options,"partitions") : DDLFamilyProfiles.POSITIONAL_QUADRANTS;
		final int overlap=options.containsKey("overlap") ? integer(options,"overlap") : DDLFamilyProfiles.POSITIONAL_OVERLAP;
		if(positionMode.positional()){DDLFamilyProfiles.validatePositionalLayout(partitions,overlap,geometry,positionMode);}
		require(positionMode.positional() || !(options.containsKey("partitions") || options.containsKey("overlap")),
			"partitions/overlap require positional reference construction");
		final long expectedAssigned=number(options,"expectedmembers"), maxMembers=number(options,"maxbucketmembers");
		final long maxBytes=number(options,"maxbucketbytes"), maxOutput=number(options,"maxoutputbytes");
		require(nf>0 && batches>0 && buckets>0 && bucket<buckets && maxReference>0 && maxReference<=100000
			&& (long)batches*buckets<=Integer.MAX_VALUE && maxMembers>0 && maxBytes>0 && maxOutput>0,
			"Invalid dimensions or resource caps");
		final Path out=Paths.get(required(options,"out")).toAbsolutePath().normalize();
		require(!Files.exists(out) && Files.isDirectory(out.getParent()),"Output must be fresh under an existing directory");
		final Path families=input(options,"families"), chunks=input(options,"chunks");
		require(options.containsKey("models")==options.containsKey("modelssha80"),"models and modelssha80 are required together");
		require(options.containsKey("consensus")==options.containsKey("consensussha80"),"consensus and consensussha80 are required together");
		final boolean sharedConsensus=options.containsKey("consensus");
		require(!references || !sharedConsensus, "consensus= applies only to lower calibration");
		require(references || sharedConsensus!=options.containsKey("models"), "Lower calibration requires exactly one of models= or consensus=");
		final Path models=references || sharedConsensus ? null : input(options,"models");
		final Path consensus=sharedConsensus ? input(options,"consensus") : null;
		Shared.AMINO_IN=true; Shared.setThreads(threads); ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final DDLPentamerBlacklist blacklist=options.containsKey("blacklist") ?
			DDLPentamerBlacklist.load(input(options,"blacklist"),required(options,"blacklistsha80")) : null;
		final long started=System.nanoTime();
		final Family[] rows=loadFamilies(families,nf,expectedAssigned);
		if(!options.containsKey("familyranks")){options.put("familyranks","all");}
		selectFamilies(options.get("familyranks"),rows,bucket,buckets);
		long selected=0;
		for(Family f : rows){if(f.selected){selected+=f.assigned;}}
		require(selected<=maxMembers && Math.multiplyExact(selected,24L)<=maxBytes,
			"Bucket exceeds member or score-array/reference-array byte cap");
		final ArrayList<Chunk> parts=loadChunks(chunks,rows,bucket,buckets,batches,expectedAssigned,memberRange);
		if(references){buildReferences(options,rows,parts,bucket,buckets,maxBytes,maxOutput,out,started,threads,sketch,partitions,overlap,blacklist,geometry); return;}
		if(consensus==null){loadModels(models,rows,bucket,buckets,maxReference);}
		else{loadConsensus(consensus,rows,maxReference);}
		loadReferences(options,rows,bucket,buckets);
		readMembers(parts,rows,bucket,buckets,maxBytes,true,retainScores,memberRange);
		final long loadEnd=System.nanoTime();
		for(Family f : rows){if(f.spectrum!=null){f.spectrum.validate(); f.norm=ProteinSpectrum.norm(f.spectrum);}}
		final long spectrumEnd=System.nanoTime();
		if(retainScores){Files.createDirectory(out);}
		final int scoreShards=runFamilies(rows,threads,false,maxOutput,retainScores ? out : null);
		final long scoreEnd=System.nanoTime();
		for(String key : new String[]{"families","chunks","references"}){input(options,key);}
		if(consensus==null){input(options,"models");}else{input(options,"consensus");}
		if(!retainScores){Files.createDirectory(out);}
		final ByteStreamWriter cut=writer(out.resolve("cutoffs.tsv")), spectra=writer(out.resolve("spectra.tsv"));
		final ByteStreamWriter oldCut=retainScores ? writer(out.resolve("cutoffs99.tsv")) : null;
		final ByteBuilder line=new ByteBuilder(); int count=0,empty=0;
		try{
			cut.println(CUTOFF_HEADER); spectra.println(ProteinSpectrum.header(2));
			if(oldCut!=null){oldCut.println(CUTOFF_HEADER);}
			for(Family f : rows){
				if(!f.selected){continue;}
				line.clear().append(f.rank).tab().append(f.id).tab().append(f.assigned).tab();
				if(f.assigned==0){line.append("NO_MEMBERS\tNA\t0\tNA\t0\tNA\t0"); empty++;}
				else{
					require(f.scored==f.assigned,"Missing assigned-member scores: "+f.id);
					line.append("CALIBRATED");
					for(Quantile q : f.cutoffs){
						line.tab().append(Double.toString(q.cutoff)).tab().append(q.retained);
					}
				}
				cut.println(line);
				if(oldCut!=null){
					if(f.assigned==0){oldCut.println(line);}
					else{
						line.clear().append(f.rank).tab().append(f.id).tab().append(f.assigned).append("\tCALIBRATED");
						for(Quantile q:f.oldCutoffs){line.tab().append(Double.toString(q.cutoff)).tab().append(q.retained);}
						oldCut.println(line);
					}
				}
				if(f.spectrum!=null){ProteinSpectrum.writeRow(line,"tail",f.id,f.spectrum); spectra.print(line);}
				count++;
			}
		}finally{require(!(cut.poisonAndWait() | spectra.poisonAndWait() | (oldCut!=null && oldCut.poisonAndWait())),"Calibration output writer failed");}
		final long written=System.nanoTime();
		writeTimings(out,rows);
		write(out.resolve("summary.tsv"),"bucket\tbuckets\tfamilies\tzero_families\tmembers\talignments\tload_seconds\tspectrum_seconds\tscore_quantile_seconds\twrite_seconds\n"+
			bucket+"\t"+buckets+"\t"+count+"\t"+empty+"\t"+selected+"\t"+selected+"\t"+
			seconds(started,loadEnd)+"\t"+seconds(loadEnd,spectrumEnd)+"\t"+seconds(spectrumEnd,scoreEnd)+"\t"+seconds(scoreEnd,written)+"\n");
		write(out.resolve("recipe.txt"),String.join("\n",args)+"\nschema=family_all_member_calibration_v2\n"+
			(retainScores ? "retention=cosine99_identity99.75_blosum99.5\ncutoff=sorted_ascending_index_floor_N_over_100_400_200;small_families_use_minimum\nraw_scores=original_double_identity_and_blosum;per_worker_shards;family_and_member_ID_retained\n" :
			"retention=99/100\ncutoff=sorted_ascending_index_floor_members_over_100\n")+"inclusive=true\n"+
			"population=all_winner_assigned_members_no_holdout_selection\nreference=fixed_truth_label_consensus\ncosine=all_member_dimer_counts_no_added_consensus\n"+
			"identity=100_matches_over_all_alignment_columns_including_gaps\nblosum=A_over_full_query_self_score\nalignment=GlocalAminoLinear_BLOSUM62_linear_gap4_query_global_reference_local\n"+
			"family_workers="+threads+"\n");
		final ArrayList<String> outputFiles=new ArrayList<String>(Arrays.asList("cutoffs.tsv","spectra.tsv","summary.tsv","recipe.txt","timing.tsv"));
		if(retainScores){
			outputFiles.add("cutoffs99.tsv");
			for(int i=0;i<scoreShards;i++){outputFiles.add("raw_scores_"+i+".tsv");}
		}
		long bytes=0;
		for(String name:outputFiles){bytes=Math.addExact(bytes,Files.size(out.resolve(name)));}
		require(bytes<=maxOutput,"Calibration output byte cap exceeded; no PASS");
		write(out.resolve("PASS"),"FAMILY_ALL_MEMBER_CALIBRATION_PASS\n");
		outputFiles.add("PASS");
		final StringBuilder pins=new StringBuilder("file\tsha80\n");
		for(String name:outputFiles){
			pins.append(name).append('\t').append(DigestSuffix.file(out.resolve(name).toString())).append('\n');
		}
		write(out.resolve("result.sha80.tsv"),pins.toString());
		System.out.println("FAMILY_ALL_MEMBER_CALIBRATION_PASS families="+count+" members="+selected+" alignments="+selected);
	}

	/** One bounded AA load; at most t live family vote maps, discarded before each next family. */
	private static void buildReferences(HashMap<String,String> options,Family[] rows,ArrayList<Chunk> parts,int bucket,int buckets,
			long cap,long maxOutput,Path out,long started,int threads,String sketch,int partitions,int overlap,DDLPentamerBlacklist blacklist,
			ProteinKmerGeometry geometry) throws Exception{
		long members=0;
		for(Family f : rows){if(f.selected){members+=f.assigned;}}
		final MemberLengthRange memberRange=MemberLengthRange.from(options);
		readMembers(parts,rows,bucket,buckets,cap,false,false,memberRange);
		final long loaded=System.nanoTime();
		runFamilies(rows,threads,true,maxOutput,null,sketch,partitions,overlap,blacklist,geometry);
		Files.createDirectory(out);
		final ByteStreamWriter profiles=writer(out.resolve("profiles.tsv")), spectra=writer(out.resolve("spectra.tsv")), census=writer(out.resolve("members.tsv"));
		final ByteBuilder line=new ByteBuilder();
		int nonempty=0,empty=0; long residues=0,words=0;
		try{
			profiles.println(DDLFamilyProfiles.HEADER16); spectra.println(ProteinSpectrum.header(2));
			census.println("roster_rank\tfamily_id\tmembers\tresidues\tvalid_words");
			for(Family f : rows){
				if(!f.selected){continue;}
				if(f.assigned==0){empty++;}
				else{
					profiles.print(f.profileRows); f.profileRows=null;
					ProteinSpectrum.writeRow(line,"tail",f.id,f.spectrum); spectra.print(line); nonempty++;
				}
				residues+=f.bases; words+=f.words;
				line.clear().append(f.rank).tab().append(f.id).tab().append(f.assigned).tab().append(f.bases).tab().append(f.words);
				census.println(line);
			}
		}finally{require(!(profiles.poisonAndWait() | spectra.poisonAndWait() | census.poisonAndWait()),"Reference output writer failed");}
		for(String key : new String[]{"families","chunks"}){input(options,key);}
		writeTimings(out,rows);
		final ProteinPositionMode positionMode=ProteinPositionMode.fromSketch(sketch);
		final boolean positional=positionMode.positional();
		final boolean legacyK5=geometry.isExact(5);
		write(out.resolve("summary.tsv"),"namespace\tfamilies_sha80\tchunks_sha80\tk\tseed\texponent\tregisterbits\tbuckets\tshard\tshards\tmanifest_families\tbuilt_families\tzero_families\tmembers\tresidues\tvalid_words\tfamilyranks\tload_seconds\tbuild_seconds\tsketch"+
			(legacyK5 ? "" : "\tgeometry")+
			DDLFamilyProfiles.positionMetadataColumns(positionMode)+
			(blacklist==null ? "" : "\tblacklist_sha80\tblacklist_words")+"\n"+
			"tail\t"+options.get("familiessha80")+"\t"+options.get("chunkssha80")+"\t"+geometry.k()+"\t0\t6\t16\t1024\t"+bucket+"\t"+buckets+"\t"+rows.length+"\t"+nonempty+"\t"+empty+"\t"+members+"\t"+residues+"\t"+words+"\t"+options.get("familyranks")+"\t"+seconds(started,loaded)+"\t"+seconds(loaded,System.nanoTime())+"\t"+sketch+
			(legacyK5 ? "" : "\t"+geometry.name())+
			DDLFamilyProfiles.positionMetadataValues(positionMode, partitions, overlap, geometry)+
			(blacklist==null ? "" : "\t"+blacklist.sha80()+"\t"+blacklist.size())+"\n");
		final StringBuilder recipe=new StringBuilder("schema=family_full_references_v1\npopulation=all_winner_assigned_members\nprofile=")
			.append(legacyK5 ? "k5_width1024_legacy16_seed0" : geometry.name()+"_width1024_a48_seed0")
			.append("\ngeometry=").append(geometry.name()).append("\nsketch=").append(sketch)
			.append("\ncosine=all_member_dimer_counts_no_added_consensus\n");
		if(blacklist!=null){
			input(options,"blacklist");
			recipe.append("blacklist=").append(required(options,"blacklist")).append("\nblacklistsha80=")
				.append(blacklist.sha80()).append("\nblacklist_words=").append(blacklist.size()).append('\n');
		}
		if(positional){
			recipe.append("partitions=").append(partitions).append("\noverlap=").append(overlap)
				.append("\ntag_bits=").append(DDLFamilyProfiles.positionalTagBits(partitions));
			if(positionMode.fixedCoordinates()){recipe.append("\nsegment_words=").append(DDLFamilyProfiles.positionalSegmentWords(partitions,geometry)).append('\n');}
			else{recipe.append("\nposition_basis=residue-length\n");}
		}
		for(String key : new String[]{"mode","families","familiessha80","chunks","chunkssha80","bucket","buckets","familyranks","expectedfamilies","expectedmembers","expectedbatches","maxbucketmembers","maxbucketbytes","maxoutputbytes","t","out"}){
			recipe.append(key).append('=').append(required(options,key)).append('\n');
		}
		recipe.append(memberRange.recipe("="));
		write(out.resolve("recipe.txt"),recipe.toString());
		long bytes=0; for(String name : REFERENCE_FILES){if(!name.equals("PASS")){bytes=Math.addExact(bytes,Files.size(out.resolve(name)));}}
		require(bytes<=maxOutput,"Reference output cap exceeded; no PASS");
		write(out.resolve("PASS"),"FAMILY_FULL_REFERENCES_PASS\n");
		final StringBuilder pins=new StringBuilder("file\tsha80\n");
		for(String name : REFERENCE_FILES){pins.append(name).append('\t').append(DigestSuffix.file(out.resolve(name).toString())).append('\n');}
		write(out.resolve("result.sha80.tsv"),pins.toString());
		System.out.println("FAMILY_FULL_REFERENCES_PASS members="+members+" profiles="+nonempty+" zero_families="+empty);
	}

	private static void appendProfile(ByteBuilder line,Family family,DDLFamilyProfiles.Profile profile,long bases,long words){
		require(profile.members()==family.assigned && profile.buckets()==1024,"Profile metadata must describe every assigned member");
		for(boolean literal : new boolean[]{true,false}){
			line.append("tail\t").append(family.id).append("\t1024\t").append(family.assigned).tab().append(bases).tab().append(words).tab().append(literal ? "literal" : "nonempty");
			for(int kind=0;kind<3;kind++){
				line.tab(); for(int b=0;b<1024;b++){if(b>0){line.append(',');} line.append(kind==0 ? (int)profile.value(b,literal) : kind==1 ? profile.votes(b,literal) : profile.occupied(b));}
			}
			line.nl();
		}
	}

	private static void loadReferences(HashMap<String,String> options,Family[] rows,int bucket,int buckets) throws Exception{
		final Path manifest=input(options,"references"), dir=manifest.getParent();
		final java.util.List<String> lines=Files.readAllLines(manifest,StandardCharsets.US_ASCII);
		require(lines.size()==REFERENCE_FILES.length+1 && lines.get(0).equals("file\tsha80"),"Require exact full-reference manifest");
		final HashMap<String,Path> files=new HashMap<String,Path>();
		for(int i=1;i<lines.size();i++){
			final String[] fields=lines.get(i).split("\t",-1);
			require(fields.length==2,"Reference manifest field count");
			final String logical=logicalReferenceName(fields[0]);
			require(files.put(logical,dir.resolve(fields[0]))==null,"Duplicate reference manifest logical file: "+logical);
			pin(dir.resolve(fields[0]),fields[1]);
		}
		for(String name : REFERENCE_FILES){require(files.containsKey(name),"Missing reference manifest file: "+name);}
		require(new String(Files.readAllBytes(files.get("PASS")),StandardCharsets.US_ASCII).equals("FAMILY_FULL_REFERENCES_PASS\n"),"Full-member references have not passed");
		final HashMap<String,String> recipe=new HashMap<String,String>();
		for(String line : Files.readAllLines(files.get("recipe.txt"),StandardCharsets.US_ASCII)){
			final String[] pair=Parse.splitOnFirst(line,'='); require(pair.length==2 && recipe.put(pair[0],pair[1])==null,"Malformed reference recipe");
		}
		final MemberLengthRange memberRange=MemberLengthRange.from(options);
		require(memberRange.same(MemberLengthRange.from(recipe)), "Reference member range differs from calibration inputs");
		require("family_full_references_v1".equals(recipe.get("schema")) && "all_winner_assigned_members".equals(recipe.get("population"))
			&& "all_member_dimer_counts_no_added_consensus".equals(recipe.get("cosine")),"Wrong reference population/kernel recipe");
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(options.containsKey("k") ? integer(options,"k") : 5,
			options.containsKey("aabits") ? options.get("aabits") : options.get("aapartition"));
		final boolean legacyK5=geometry.isExact(5);
		require((legacyK5 ? "k5_width1024_legacy16_seed0" : geometry.name()+"_width1024_a48_seed0").equals(recipe.get("profile")) &&
			geometry.name().equals(recipe.containsKey("geometry") ? recipe.get("geometry") : ProteinKmerGeometry.exact(5).name()),
			"Wrong reference protein k-mer geometry");
		final String sketch=recipe.containsKey("sketch") ? recipe.get("sketch") : DDLFamilyProfiles.SKETCH_PLAIN;
		require(DDLFamilyProfiles.SKETCH_PLAIN.equals(sketch),"Lower scoring requires plain DDL references, not "+sketch);
		for(String key : new String[]{"familiessha80","chunkssha80","bucket","buckets","familyranks","expectedfamilies","expectedmembers","expectedbatches"}){
			require(required(options,key).equals(recipe.get(key)),"Reference input binding differs: "+key);
		}
		final LinkedHashMap<String,ProteinSpectrum.Spectrum> spectra=ProteinSpectrum.load(files.get("spectra.tsv").toString());
		final java.util.List<String> census=Files.readAllLines(files.get("members.tsv"),StandardCharsets.US_ASCII); int index=1;
		require(!census.isEmpty() && census.get(0).equals("roster_rank\tfamily_id\tmembers\tresidues\tvalid_words"),"Wrong reference census header");
		for(Family f : rows){
			if(!f.selected){continue;}
			require(index<census.size(),"Missing reference census family");
			final String[] fields=census.get(index++).split("\t",-1);
			require(fields.length==5 && Integer.parseInt(fields[0])==f.rank && fields[1].equals(f.id) && Integer.parseInt(fields[2])==f.assigned,"Reference census differs from assigned roster");
			final long bases=Long.parseLong(fields[3]), words=Long.parseLong(fields[4]);
			require(memberRange.containsTotal(f.assigned,bases) && words>=0 &&
				words<=bases-(geometry.span()-1L)*f.assigned,"Reference census residues/words invalid");
			if(f.assigned==0){continue;}
			f.spectrum=spectra.remove("tail:"+f.id);
			require(f.spectrum!=null && f.spectrum.k==2 && f.spectrum.total<=bases-f.assigned,"Missing/invalid full-member dimer spectrum: "+f.id);
		}
		require(spectra.isEmpty() && index==census.size(),"Extra family spectra/census rows outside bound population");
		input(options,"references");
	}

	private static String logicalReferenceName(String name){
		if(name.equals("profiles.tsv.gz")){return "profiles.tsv";}
		if(name.equals("spectra.tsv.gz")){return "spectra.tsv";}
		for(String allowed : REFERENCE_FILES){if(name.equals(allowed)){return name;}}
		throw new IllegalArgumentException("Unexpected reference manifest file: "+name);
	}

	private static void selectFamilies(String selection,Family[] rows,int bucket,int buckets){
		require(selection!=null && !selection.isEmpty(),"Empty family selection");
		if(selection.equals("all")){for(Family f : rows){f.selected=f.rank%buckets==bucket;} return;}
		int previous=-1;
		for(String token : selection.split(",",-1)){
			final int rank=Integer.parseInt(token);
			require(rank>previous && rank<rows.length && rank>=0 && rank%buckets==bucket && token.equals(Integer.toString(rank)),"Explicit family ranks must be canonical ascending distinct ranks in the selected bucket");
			rows[rank].selected=true; previous=rank;
		}
	}

	private static void writeTimings(Path out,Family[] rows) throws Exception{
		final StringBuilder table=new StringBuilder("roster_rank\tfamily_id\tmembers\treference_build_wall_seconds\tpositive_score_wall_seconds\tquantile_wall_seconds\treference_build_cpu_seconds\tpositive_score_cpu_seconds\tquantile_cpu_seconds\n");
		for(Family f : rows){if(f.selected){table.append(f.rank).append('\t').append(f.id).append('\t').append(f.assigned).append('\t')
			.append(f.buildNanos*1e-9).append('\t').append(f.scoreNanos*1e-9).append('\t').append(f.quantileNanos*1e-9).append('\t')
			.append(cpuSeconds(f.buildCpuNanos)).append('\t').append(cpuSeconds(f.scoreCpuNanos)).append('\t').append(cpuSeconds(f.quantileCpuNanos)).append('\n');}}
		write(out.resolve("timing.tsv"),table.toString());
	}

	/** Java8-safe current-worker CPU, not process CPU; never enable a VM-wide counter implicitly. */
	public static long cpuTime(ThreadMXBean bean){
		try{return bean!=null && bean.isCurrentThreadCpuTimeSupported() && bean.isThreadCpuTimeEnabled() ? bean.getCurrentThreadCpuTime() : -1;}
		catch(SecurityException | UnsupportedOperationException unavailable){return -1;}
	}
	/** Unavailable or discontinuous counters cannot establish a compute baseline. */
	public static long cpuElapsed(long start,long end){return start<0 || end<start ? -1 : end-start;}
	private static String cpuSeconds(long nanos){return nanos<0 ? "NA" : Double.toString(nanos*1e-9);}
	private static ThreadMXBean cpuBean(){
		try{return ManagementFactory.getThreadMXBean();}catch(SecurityException | UnsupportedOperationException unavailable){return null;}
	}

	/** Member retention is independent of whether a caller writes raw scores or comparison cutoffs. */
	public enum MemberRetention{
		/** Retain at least 99% for all three metrics. */
		UNIFORM_99,
		/** Retain at least 99% cosine, 99.75% identity and 99.5% BLOSUM. */
		METRIC_SPECIFIC
	}

	/** Fits cosine, identity and BLOSUM cutoffs in that order, sorting each score row in place once.
	 * Rows must contain the same nonzero member count and finite scores. The returned array is caller-owned.
	 * If nonnull, legacyCutoffs must have three slots; they receive the uniform 99% comparison from the same sort.
	 * This method performs no I/O and retains no score arrays. A failed call may leave earlier rows sorted.
	 */
	public static Quantile[] memberCutoffs(double[][] scores,MemberRetention retention,Quantile[] legacyCutoffs){
		require(scores!=null && scores.length==3 && retention!=null,"Member cutoffs require three metric rows and an explicit retention policy");
		require(scores[0]!=null && scores[0].length>0,"A zero family has no fitted member cutoffs");
		for(double[] row:scores){require(row!=null && row.length==scores[0].length,"All three metrics must describe the same member population");}
		require(legacyCutoffs==null || legacyCutoffs.length==3,"Legacy comparison output must have one slot per metric");
		final Quantile[] cutoffs=new Quantile[3];
		for(int metric=0;metric<3;metric++){
			final int divisor=retentionDivisor(retention,metric);
			cutoffs[metric]=quantile(scores[metric],divisor);
			if(legacyCutoffs!=null){legacyCutoffs[metric]=divisor==100 ? cutoffs[metric] : sortedQuantile(scores[metric],100);}
		}
		return cutoffs;
	}

	/** Fits per-metric cutoffs when alignment metrics are a measured prefix, not every assigned member.
	 * Failed alignments consume the same rejected-member quota that low measured scores would have used.
	 * The supplied rows are sorted only over their measured prefixes; tail slots remain caller scratch.
	 */
	public static MemberCutoff[] measuredMemberCutoffs(double[][] scores,int[] measuredCounts,
			int assigned,MemberRetention retention){
		require(scores!=null && scores.length==3 && measuredCounts!=null && measuredCounts.length==3 && retention!=null,
			"Measured cutoffs require three metric rows, three counts and an explicit retention policy");
		final MemberCutoff[] cutoffs=new MemberCutoff[3];
		for(int metric=0;metric<3;metric++){
			cutoffs[metric]=measuredMemberCutoff(scores[metric],measuredCounts[metric],assigned,
				retentionDivisor(retention,metric));
		}
		return cutoffs;
	}

	/** One lower cutoff from finite prefix values plus the full assigned-member denominator. */
	public static MemberCutoff measuredMemberCutoff(double[] scores,int measured,int assigned,int divisor){
		require(assigned>0 && divisor>1,"Measured cutoff requires assigned members and a retained cutoff divisor");
		require(measured>=0 && measured<=assigned,"Measured count must be within the assigned member count");
		require(scores!=null && measured<=scores.length,"Measured prefix must fit in the supplied score row");
		for(int i=0;i<measured;i++){require(Double.isFinite(scores[i]),"Nonfinite measured metric cannot become a fitted cutoff");}
		final int failed=assigned-measured, quota=assigned/divisor;
		if(measured==0){return MemberCutoff.unavailable(assigned,failed,quota);}
		Arrays.sort(scores,0,measured);
		final int measuredRejections=Math.max(0,quota-failed);
		final double cutoff=scores[measuredRejections];
		int rejected=0;
		while(rejected<measured && scores[rejected]<cutoff){rejected++;}
		return MemberCutoff.measured(cutoff,assigned,failed,quota,measured-rejected,failed<=quota);
	}

	private static int retentionDivisor(MemberRetention retention,int metric){
		return retention==MemberRetention.UNIFORM_99 || metric==0 ? 100 : metric==1 ? 400 : 200;
	}

	/** Exact lower empirical quantile; ties can only increase retained count. Mutates the supplied array. */
	public static Quantile quantile(double[] scores){
		return quantile(scores,100);
	}
	/** Inclusive empirical lower; N<divisor uses the minimum, including the 400/200 member boundaries. */
	public static Quantile quantile(double[] scores,int divisor){
		require(scores!=null && scores.length>0,"A zero family has no fitted cutoff");
		require(divisor>1,"Quantile divisor must exceed one so a retained cutoff exists");
		for(double value : scores){require(Double.isFinite(value),"Nonfinite member metric cannot be silently discarded");}
		Arrays.sort(scores);
		return sortedQuantile(scores,divisor);
	}
	private static Quantile sortedQuantile(double[] scores,int divisor){
		assert(scores.length>0 && divisor>0) : "Quantile input was validated before sorting";
		final double cutoff=scores[scores.length/divisor];
		int rejected=0;
		while(rejected<scores.length && scores[rejected]<cutoff){rejected++;}
		final int retained=scores.length-rejected;
		require(retained>=scores.length-scores.length/divisor,"Inclusive quantile must retain the exact integer quota");
		return new Quantile(cutoff,retained);
	}

	public static final class Quantile{
		private Quantile(double value,int count){cutoff=value; retained=count;}
		public final double cutoff; public final int retained;
	}

	public static final class MemberCutoff{
		private MemberCutoff(State state_,double cutoff_,int assigned_,int failed_,int quota_,int retained_,boolean achievable_){
			require(state_!=null && Double.isFinite(cutoff_) && assigned_>0 && failed_>=0 &&
				failed_<=assigned_ && quota_>=0 && retained_>=0 && retained_<=assigned_-failed_,
				"Member cutoff accounting must be finite and internally bounded");
			require(state_==State.MEASURED ? retained_>0 && failed_<assigned_ :
				retained_==0 && failed_==assigned_ && !achievable_,
				"Unavailable member cutoffs must be explicit and measured cutoffs need measured values");
			state=state_; cutoff=cutoff_; assigned=assigned_; failed=failed_;
			measured=assigned-failed; rejectionQuota=quota_; retained=retained_; quotaAchievable=achievable_;
		}
		private static MemberCutoff measured(double cutoff,int assigned,int failed,int quota,int retained,boolean achievable){
			return new MemberCutoff(State.MEASURED,cutoff,assigned,failed,quota,retained,achievable);
		}
		private static MemberCutoff unavailable(int assigned,int failed,int quota){
			return new MemberCutoff(State.UNAVAILABLE,0,assigned,failed,quota,0,false);
		}
		public boolean measured(){return state==State.MEASURED;}
		public double cutoff(){
			if(!measured()){throw new IllegalArgumentException("No measured values are available for this member cutoff");}
			return cutoff;
		}
		public enum State{MEASURED,UNAVAILABLE}
		public final State state;
		public final int assigned,measured,failed,rejectionQuota,retained;
		public final boolean quotaAchievable;
		private final double cutoff;
	}

	/** Shared metric definition for fixtures/consumer. Output is caller-owned; the aligner returns one alignment. */
	public static void measure(ProteinSequence query,ProteinSpectrum.Spectrum querySpectrum,ProteinSequence reference,
			ProteinSpectrum.Spectrum familySpectrum,double familyNorm,Metrics out){
		require(query!=null && querySpectrum!=null,"Query inputs required");
		measure(query,querySpectrum,reference,familySpectrum,ProteinSpectrum.norm(querySpectrum),familyNorm,ProteinSearcher.selfScore(query),out);
	}

	/** Multi-candidate callers compute query norm and full-query self score once, family norms once per library. */
	public static void measure(ProteinSequence query,ProteinSpectrum.Spectrum querySpectrum,ProteinSequence reference,
			ProteinSpectrum.Spectrum familySpectrum,double queryNorm,double familyNorm,int querySelf,Metrics out){
		require(querySpectrum!=null && familySpectrum!=null && out!=null && querySpectrum.k==2 && familySpectrum.k==2,"Dimer spectra and metric output required");
		out.cosine=ProteinSpectrum.cosine(querySpectrum,familySpectrum,queryNorm,familyNorm);
		require(Double.isFinite(out.cosine) && out.cosine>=0 && out.cosine<=1,"Invalid dimer cosine");
		align(query,reference,querySelf,out);
	}

	/** For survivors of an already computed cosine gate; preserves out.cosine and never recomputes its dot product. */
	public static void align(ProteinSequence query,ProteinSequence reference,int querySelf,Metrics out){
		//TODO: Caller contract mismatch: this charges query prefixes and leaves the reference end free,
		// unlike the requested left-clipped/right-pinned tail comparison (SHORT_MISS_ANATOMY_20261104).
		// Repair caller and calibration together; existing member/negative gates use this full-query metric.
		require(query!=null && reference!=null && out!=null && query.enc.length>0 && reference.enc.length>0,"Nonempty alignment inputs required");
		require(querySelf>0,"A/Sq requires positive full-query self score");
		final AAAlignment alignment=GlocalAminoLinear.align(query.enc,reference.enc,false);
		out.identityPct=alignment.pident(); out.blosum=alignment.rawScore/(double)querySelf;
		require(Double.isFinite(out.identityPct)
			&& out.identityPct>=0 && out.identityPct<=100 && Double.isFinite(out.blosum),"Invalid calibrated metric");
	}

	public static final class Metrics{public double cosine,identityPct,blosum;}

	/** Input identity policy; structural validation never establishes byte identity. */
	public enum VerificationMode{
		/** Recompute the retained sha80 before and after reading, as existing callers require. */
		SHA80,
		/** Validate structure and metadata stability; declared identity may be inherited or unmeasured. */
		STRUCTURAL
	}

	/** Preparation inputs sharing the canonical CLI readers. One coordinator loads each
	 * bucket once; returned member/consensus arrays are detached and owned by the caller.
	 * The caller releases a bucket after preparation and calls verifyComplete before publishing.
	 * Decoded/retained member caps exclude other workers, prepared features, panels and JVM overhead.
	 */
	public static final class MemberInputs{
		public MemberInputs(Path families, String familyPin, Path chunks, String chunkPin,
				Path consensus, String consensusPin, int expectedFamilies, long expectedMembers,
				int buckets_, int batches_, int maxReference){
			this(families, familyPin, chunks, chunkPin, consensus, consensusPin, expectedFamilies, expectedMembers,
				buckets_, batches_, maxReference, MemberLengthRange.LEGACY);
		}
		public MemberInputs(Path families, String familyPin, Path chunks, String chunkPin,
				Path consensus, String consensusPin, int expectedFamilies, long expectedMembers,
				int buckets_, int batches_, int maxReference, MemberLengthRange memberRange_){
			this(families, familyPin, chunks, chunkPin, consensus, consensusPin, expectedFamilies, expectedMembers,
				buckets_, batches_, maxReference, memberRange_, VerificationMode.SHA80);
		}
		/** Explicit structural mode preserves reader invariants without measuring byte identity.
		 * Callers must record that limitation and supply their independent population control.
		 * Existing constructors retain SHA80 verification; this policy is local to this session.
		 */
		public MemberInputs(Path families, String familyPin, Path chunks, String chunkPin,
				Path consensus, String consensusPin, int expectedFamilies, long expectedMembers,
				int buckets_, int batches_, int maxReference, MemberLengthRange memberRange_, VerificationMode verificationMode_){
			require(memberRange_!=null, "Preparation inputs require their member-length contract");
			require(verificationMode_!=null, "Preparation inputs require an explicit input verification policy");
			memberRange=memberRange_;
			verificationMode=verificationMode_;
			structural=verificationMode==VerificationMode.STRUCTURAL ? new StructuralInputs() : null;
			require(families!=null && chunks!=null && consensus!=null && familyPin!=null && chunkPin!=null && consensusPin!=null,
				"Preparation requires the family, chunk and consensus paths with their fingerprints");
			require(Shared.AMINO_IN, "Configure amino-acid input before opening preparation inputs");
			require(expectedFamilies>0 && expectedMembers>=0 && buckets_>0 && batches_>0 &&
				(long)buckets_*batches_<=Integer.MAX_VALUE && maxReference>0 && maxReference<=1000,
				"Preparation dimensions must fit the canonical chunk grid and1..1000-AA consensus contract");
			paths=new Path[]{families, chunks, consensus}; pins=new String[]{familyPin, chunkPin, consensusPin};
			for(int i=0; i<paths.length; i++){verifyInput(paths[i], pins[i], verificationMode, structural);}
			buckets=buckets_; batches=batches_; members=expectedMembers; loaded=new boolean[buckets];
			rows=loadFamilies(families, expectedFamilies, expectedMembers);
			for(Family row:rows){row.selected=true;}
			loadConsensus(consensus, rows, maxReference, true);
			// Preparation encodes the original bytes itself; do not retain a redundant encoded roster.
			for(Family row:rows){row.reference=null;}
			for(int i=0; i<paths.length; i++){verifyInput(paths[i], pins[i], verificationMode, structural);}
		}
		public int familyCount(){return rows.length;}
		/** Preserves ascending physical ranks and canonical per-family order. A failed read invalidates
		 * this input session so a partial bucket cannot later be mistaken for completed preparation.
		 */
		public MemberInput[] loadBucket(int bucket, long maxMembers, long maxBytes){
			require(!failed && bucket>=0 && bucket<buckets && !loaded[bucket], "Preparation bucket must be in range, unread and from a valid input session");
			require(maxMembers>0 && maxBytes>0, "Preparation member and byte limits must be positive");
			try{
				long selected=0; int count=0;
				for(Family row:rows){
					row.selected=row.rank%buckets==bucket;
					if(row.selected){selected=Math.addExact(selected, row.assigned); count++;}
				}
				require(selected<=maxMembers && Math.multiplyExact(selected, 24L)<=maxBytes,
					"Preparation bucket exceeds the member or reserved-score byte cap");
				verifyInput(paths[1], pins[1], verificationMode, structural);
				final ArrayList<Chunk> parts=loadChunks(paths[1], rows, bucket, buckets, batches, members, memberRange, verificationMode);
				decodedMemberBytes=Math.addExact(decodedMemberBytes,
					readMembers(parts, rows, bucket, buckets, maxBytes, true, false, memberRange, verificationMode, structural));
				verifyInput(paths[1], pins[1], verificationMode, structural);
				final MemberInput[] result=new MemberInput[count]; int slot=0;
				for(Family row:rows){
					if(!row.selected){continue;}
					result[slot++]=new MemberInput(row.rank, row.id, row.sequences, row.consensusRaw);
					row.sequences=null; row.consensusRaw=null;
				}
				assert(slot==count) : "Every selected family, including empty families, must transfer exactly once";
				loaded[bucket]=true;
				return result;
			}catch(RuntimeException | Error failure){failed=true; throw failure;}
		}
		public void verifyComplete(){
			require(!failed, "Cannot publish preparation after a failed input read");
			for(int bucket=0; bucket<buckets; bucket++){require(loaded[bucket], "Preparation has not consumed bucket"+bucket);}
			verifyIdentities();
		}
		/** Verify exactly the declared bucket set for a cost pilot, without claiming full payload coverage.
		 * Original family ranks and full roster dimensions remain unchanged; all read inputs are bracketed.
		 */
		public void verifySelectedBuckets(int[] selected){
			require(!failed, "Cannot publish a selected-bucket pilot after a failed input read");
			require(selected!=null && selected.length>0, "A pilot must declare a nonempty bucket set");
			int previous=-1;
			for(int bucket:selected){
				require(bucket>previous && bucket<buckets, "Pilot buckets must be ascending, distinct and in the original grid");
				previous=bucket;
			}
			int next=0;
			for(int bucket=0; bucket<buckets; bucket++){
				final boolean expected=next<selected.length && selected[next]==bucket;
				require(loaded[bucket]==expected, "Loaded bucket set differs from the declared pilot at bucket"+bucket);
				if(expected){next++;}
			}
			try{verifyIdentities();}
			catch(RuntimeException | Error failure){failed=true; throw failure;}
		}
		/** Decoded member data-record bytes only, excluding headers and line terminators; no second read. */
		public long decodedMemberRecordBytes(){
			require(!failed, "Failed input sessions cannot supply a successful decoded-byte measurement");
			return decodedMemberBytes;
		}
		/** Full and selected completion share the same identity and metadata-stability checks. */
		private void verifyIdentities(){
			for(int i=0; i<paths.length; i++){verifyInput(paths[i], pins[i], verificationMode, structural);}
			if(structural!=null){
				try{structural.verify();}
				catch(java.io.IOException | IllegalStateException failure){failed=true;throw new IllegalArgumentException("Structural input changed before completion",failure);}
			}
		}
		public final int buckets, batches;
		public final long members;
		public final MemberLengthRange memberRange;
		public final VerificationMode verificationMode;
		private final Family[] rows;
		private final Path[] paths;
		private final String[] pins;
		private final boolean[] loaded;
		private final StructuralInputs structural;
		private long decodedMemberBytes;
		private boolean failed;
	}

	/** One canonical family input, including explicit empty members and the original consensus bytes. */
	public static final class MemberInput{
		private MemberInput(int rank_, String id_, byte[][] members_, byte[] consensus_){
			require(members_!=null && consensus_!=null, "Every prepared input needs its member array and original consensus");
			rank=rank_; id=id_; members=members_; consensus=consensus_;
		}
		public final int rank;
		public final String id;
		public final byte[][] members;
		public final byte[] consensus;
	}

	private static Family[] loadFamilies(Path file,int expected,long assigned){
		final Family[] rows=new Family[expected]; final HashSet<String> ids=new HashSet<String>();
		final ByteFile input=ByteFile.makeByteFile1(file.toString(),false); final LineParser1 lp=new LineParser1('\t');
		int rank=0; long total=0;
		try{
			header(input,DDLFamilyHoldout.FAMILY_HEADER);
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line); require(lp.terms()==3 && rank<expected && lp.parseInt(0)==rank,"Noncanonical family rank/header");
				final String id=lp.parseString(1); final int count=lp.parseInt(2);
				require(id.matches("[A-Za-z0-9_.-]+") && ids.add(id) && count>=0,"Unsafe/duplicate family or negative assigned count");
				rows[rank]=new Family(rank,id,count); total=Math.addExact(total,count); rank++;
			}
		}finally{require(!input.close(),"Family metadata read failed");}
		require(rank==expected && total==assigned,"Family count/assigned conservation failed"); return rows;
	}

	private static void loadModels(Path file,Family[] rows,int bucket,int buckets,int cap){
		final ByteFile input=ByteFile.makeByteFile1(file.toString(),false); final LineParser1 lp=new LineParser1('\t'); int rank=0;
		try{
			header(input,HbmTailLibrary.MODELS_HEADER);
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line); require(lp.terms()==9 && rank<rows.length && lp.parseInt(0)==rank && lp.termEquals(rows[rank].id,1),"Fixed model roster differs from assigned family roster");
				final Family f=rows[rank++];
				if(!f.selected || f.assigned==0){continue;}
				final Path path=Paths.get(lp.parseString(2)).resolve("consensus.faa"); final String hash=lp.parseString(6);
				pin(path,hash); f.reference=consensus(path,cap,f.id); pin(path,hash);
			}
		}finally{require(!input.close(),"Fixed models read failed");}
		require(rank==rows.length,"Missing fixed consensus model rows");
	}

	private static void loadConsensus(Path file,Family[] rows,int cap){
		loadConsensus(file, rows, cap, false);
	}
	private static void loadConsensus(Path file,Family[] rows,int cap,boolean retainRaw){
		final HashMap<String,Family> byId=new HashMap<String,Family>(rows.length*2);
		int selected=0;
		for(Family f : rows){
			require(byId.put(f.id,f)==null,"Duplicate family ID in the accepted roster: "+f.id);
			if(f.selected && f.assigned>0){selected++;}
		}
		final FileFormat format=FileFormat.testInput(file.toString(),FileFormat.FASTA,null,false,true);
		final Streamer input=StreamerFactory.makeStreamer(format,null,true,-1,false,true,1);
		input.start();
		final HashSet<String> seen=new HashSet<String>();
		try{
			for(ListNum<Read> batch=input.nextList();batch!=null;batch=input.nextList()){
				for(Read read : batch.list){
					require(read!=null && read.bases!=null && read.mate==null && read.bases.length>0 && read.bases.length<=cap,
						"Shared consensus must contain bounded unpaired proteins");
					final String id=firstToken(read.id);
					final Family f=byId.get(id);
					require(f!=null && seen.add(id),"Shared consensus contains a duplicate or foreign family: "+id);
					if(retainRaw && f.selected){f.consensusRaw=Arrays.copyOf(read.bases, read.bases.length);}
					if(f.selected && f.assigned>0){f.reference=new ProteinSequence(id,read.bases); selected--;}
				}
			}
		}finally{require(!ReadWrite.closeStream(input),"Shared consensus read failed");}
		require(selected==0 && seen.size()==rows.length,"Shared consensus does not match the accepted family roster");
	}

	private static ProteinSequence consensus(Path path,int cap,String family){
		final FileFormat format=FileFormat.testInput(path.toString(),FileFormat.FASTA,null,false,true);
		final Streamer input=StreamerFactory.makeStreamer(format,null,true,-1,false,true,1); input.start();
		ProteinSequence found=null;
		try{
			for(ListNum<Read> batch=input.nextList();batch!=null;batch=input.nextList()){
				for(Read read : batch.list){
					require(read!=null && read.bases!=null && read.mate==null && read.bases.length>0 && read.bases.length<=cap && found==null,"Fixed consensus must contain exactly one bounded unpaired protein");
					// HbmTailLibrary.exportJob requires the singleton consensus ID to equal models.family_id.
					require(read.id!=null,"Missing fixed consensus header: "+path);
					int end=0; while(end<read.id.length() && !Character.isWhitespace(read.id.charAt(end))){end++;}
					final String id=read.id.substring(0,end);
					require(id.equals(family),"Fixed consensus ID "+id+" differs from models family "+family+": "+path);
					found=new ProteinSequence(id,read.bases);
				}
			}
		}finally{require(!ReadWrite.closeStream(input),"Fixed consensus read failed");}
		require(found!=null,"Missing fixed consensus sequence"); return found;
	}

	private static ArrayList<Chunk> loadChunks(Path file,Family[] rows,int bucket,int buckets,int batches,long total,MemberLengthRange memberRange){
		return loadChunks(file,rows,bucket,buckets,batches,total,memberRange,VerificationMode.SHA80);
	}
	/** Parse declared identities under the explicit session policy; never infer mode from a label. */
	private static ArrayList<Chunk> loadChunks(Path file,Family[] rows,int bucket,int buckets,int batches,long total,
			MemberLengthRange memberRange,VerificationMode mode){
		require(mode!=null,"Chunk manifest requires an explicit verification policy");
		final ArrayList<Chunk> parts=new ArrayList<Chunk>(); final HashSet<String> paths=new HashSet<String>();
		final ByteFile input=ByteFile.makeByteFile1(file.toString(),false); final LineParser1 lp=new LineParser1('\t');
		int row=0; long n=0;
		try{
			header(input,DDLFamilyHoldout.CHUNK_HEADER);
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line); require(lp.terms()==6 && row<(long)buckets*batches,"Extra/malformed chunk row");
				final int b=row%buckets;
				require(lp.termEquals(String.format(Locale.ROOT,"b%02d",row/buckets),0) && lp.parseInt(1)==b,"Chunks must preserve canonical batch/bucket order");
				final String path=lp.parseString(2),hash=lp.parseString(3); final long count=lp.parseLong(4),residues=lp.parseLong(5);
				require(Paths.get(path).isAbsolute() && paths.add(path) && validIdentity(hash,mode) && count>=0
					&& memberRange.containsTotal(count,residues),"Invalid chunk pin/path/count or member range "+memberRange);
				if(b==bucket){parts.add(new Chunk(Paths.get(path),hash,count,residues));}
				n=Math.addExact(n,count); row++;
			}
		}finally{require(!input.close(),"Chunk manifest read failed");}
		require(row==buckets*batches && n==total && parts.size()==batches,"Missing chunks or assigned count mismatch"); return parts;
	}

	/** Validate each bucket exactly once and retain canonical per-family member order.
	 * Reserve conservative AA-pointer/array and lower-score storage before starting workers.
	 * Vote maps, compact output, metadata, IO and JVM overhead are separate from this storage cap.
	 */
	private static void readMembers(ArrayList<Chunk> parts,Family[] rows,int bucket,int buckets,long cap,boolean lower,boolean retainIds,MemberLengthRange memberRange){
		readMembers(parts,rows,bucket,buckets,cap,lower,retainIds,memberRange,VerificationMode.SHA80,null);
	}
	/** Returns the already-counted decoded data-record bytes after complete bucket validation. */
	private static long readMembers(ArrayList<Chunk> parts,Family[] rows,int bucket,int buckets,long cap,boolean lower,boolean retainIds,
			MemberLengthRange memberRange,VerificationMode verificationMode,StructuralInputs structural){
		final int[] ordinals=new int[rows.length];
		final HashSet<String> genomes=new HashSet<String>(); final Buffers buffers=new Buffers();
		long bytes=0,retained=0;
		for(Family f : rows){if(f.selected){
			retained=Math.addExact(retained,24L+8L*f.assigned);
			if(lower && f.assigned>0){retained=Math.addExact(retained,120L+24L*f.assigned);}
			if(retainIds){retained=Math.addExact(retained,24L+8L*f.assigned);}
		}}
		require(retained<=cap,"AA-pointer and reserved score arrays exceed bucket cap");
		for(Family f : rows){if(f.selected){f.sequences=new byte[f.assigned][];if(retainIds){f.memberIds=new byte[f.assigned][];}}}
		for(Chunk part : parts){
			verifyInput(part.path,part.hash,verificationMode,structural);
			final ByteFile input=ByteFile.makeByteFile1(part.path.toString(),false); final LineParser1 lp=new LineParser1('\t');
			long count=0,residues=0,lastSource=-1; String genome=null;
			try{
				header(input,DDLFamilyHoldout.MEMBER_HEADER);
				for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
					lp.set(line); require(lp.terms()==6,"Joined member needs six fields");
					final int rank=lp.parseInt(0),length=lp.length(5);
					if(rank<0 || rank>=rows.length || rank%buckets!=bucket || !memberRange.contains(length)){
						throw new IllegalArgumentException("Invalid member rank/bucket/length; expected "+memberRange+", rank="+rank+" length="+length);
					}
					final Family f=rows[rank]; require(lp.termEquals(f.id,1),"Joined family identity differs from canonical roster");
					for(int field=2;field<=3;field++){
						require(lp.length(field)>0,"Empty query/genome identity");
						for(int i=0;i<lp.length(field);i++){final byte c=lp.parseByte(field,i); require(c>32 && c<127,"Noncanonical query/genome identity");}
					}
					if(genome==null || !lp.termEquals(genome,3)){
						genome=lp.parseString(3); require(genomes.add(genome),"Repeated/reordered genome in canonical chunks"); lastSource=-1;
					}
					final long source=lp.parseLong(4); require(source>=0 && source>lastSource,"Repeated/reordered source ordinal"); lastSource=source;
					require(lp.length(2)>genome.length()+1 && lp.parseByte(2,genome.length())=='|',"Query ID lacks genome namespace");
					for(int i=0;i<genome.length();i++){require(lp.parseByte(2,i)==genome.charAt(i),"Query/genome namespace differs");}
					bytes=Math.addExact(bytes,line.length); require(bytes<=cap,"Decoded member byte cap exceeded");
					require(ordinals[rank]++<f.assigned,"Too many assigned members: "+f.id);
					final byte[] raw=buffers.raw(length); System.arraycopy(line,line.length-length,raw,0,length);
					buffers.validate(raw);
					if(f.selected){
						retained=Math.addExact(retained,16L+((length+7L)&~7L)); require(retained<=cap,"Retained AA plus reserved score arrays exceed bucket cap");
						f.sequences[ordinals[rank]-1]=Arrays.copyOf(raw,length);
						if(retainIds){
							retained=Math.addExact(retained,16L+((lp.length(2)+7L)&~7L));
							require(retained<=cap,"Retained member IDs exceed bucket byte cap");
							f.memberIds[ordinals[rank]-1]=lp.parseByteArray(2);
						}
					}
					count++; residues+=length;
				}
			}finally{require(!input.close(),"Joined member read failed");}
			require(count==part.members && residues==part.residues,"Joined chunk member/residue mismatch"); verifyInput(part.path,part.hash,verificationMode,structural);
		}
		for(Family f : rows){
			if(f.rank%buckets!=bucket){continue;}
			require(ordinals[f.rank]==f.assigned,"Assigned member count not conserved: "+f.id);
		}
		return bytes;
	}

	/** Dispatch once per whole family. Workers own raw-score shards; the coordinator writes aggregate tables. */
	private static int runFamilies(Family[] rows,int threads,boolean references,long maxOutput,Path rawOut){
		return runFamilies(rows,threads,references,maxOutput,rawOut,DDLFamilyProfiles.SKETCH_PLAIN,
			DDLFamilyProfiles.POSITIONAL_QUADRANTS,DDLFamilyProfiles.POSITIONAL_OVERLAP,null);
	}

	private static int runFamilies(Family[] rows,int threads,boolean references,long maxOutput,Path rawOut,String sketchMode,int partitions,int overlap,DDLPentamerBlacklist blacklist){
		return runFamilies(rows,threads,references,maxOutput,rawOut,sketchMode,partitions,overlap,blacklist,ProteinKmerGeometry.exact(5));
	}

	private static int runFamilies(Family[] rows,int threads,boolean references,long maxOutput,Path rawOut,String sketchMode,int partitions,int overlap,DDLPentamerBlacklist blacklist,
			ProteinKmerGeometry geometry){
		require(threads>0 && maxOutput>0,"Positive worker/output bounds required");
		final ProteinPositionMode mode=ProteinPositionMode.fromSketch(sketchMode);
		if(mode.positional()){DDLFamilyProfiles.validatePositionalLayout(partitions,overlap,geometry,mode);}
		final ArrayList<Family> selected=new ArrayList<Family>();
		for(Family f : rows){if(f.selected){selected.add(f);}}
		final Family[] tasks=selected.toArray(new Family[0]);
		final AtomicInteger next=new AtomicInteger(); final AtomicLong outputBytes=new AtomicLong();
		final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>();
		for(int i=0;i<Math.min(threads,tasks.length);i++){
			workers.add(new ProcessThread(tasks,next,outputBytes,references,maxOutput,i,rawOut,sketchMode,partitions,overlap,blacklist,geometry));
		}
		try{ThreadWaiter.startThreads(workers);}
		finally{
			// A failed Thread.start must not leave earlier workers alive or join a NEW thread forever.
			final ArrayList<ProcessThread> started=new ArrayList<ProcessThread>();
			for(ProcessThread worker : workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
			ThreadWaiter.waitForThreadsToFinish(started);
		}
		for(ProcessThread worker : workers){
			if(worker.failure!=null){throw new IllegalStateException("Family worker failed: "+worker.getName(),worker.failure);}
			require(worker.success,"Family worker exited without success: "+worker.getName());
		}
		return workers.size();
	}

	/** Each family is mutated by only its claiming worker; scratch and metrics never cross workers. */
	private static final class ProcessThread extends Thread{
		ProcessThread(Family[] tasks_,AtomicInteger next_,AtomicLong bytes_,boolean references_,long maxOutput_,int id,Path rawOut,String sketchMode_,int partitions_,int overlap_,DDLPentamerBlacklist blacklist_,
				ProteinKmerGeometry geometry_){
			super("family-calibration-"+id);
			assert(tasks_!=null && next_!=null && bytes_!=null && maxOutput_>0) : "Dispatch and output bounds must outlive all family workers";
			tasks=tasks_; next=next_; outputBytes=bytes_; references=references_; maxOutput=maxOutput_;
			positionMode=ProteinPositionMode.fromSketch(sketchMode_);
			partitions=partitions_; overlap=overlap_; geometry=geometry_;
			blacklist=blacklist_;
			rawPath=rawOut==null ? null : rawOut.resolve("raw_scores_"+id+".tsv");
		}
		@Override public void run(){
			try{
				if(rawPath!=null){rawWriter=writer(rawPath);rawWriter.println(RAW_SCORE_HEADER);}
				try{
				for(int index=next.getAndIncrement();index<tasks.length;index=next.getAndIncrement()){
					final Family f=tasks[index];
					if(f.assigned>0){if(references){build(f);}else{score(f);}}
					f.sequences=null;
					f.memberIds=null;
				}
				}finally{
					if(rawWriter!=null){
						try{flushRaw();}finally{require(!rawWriter.poisonAndWait(),"Raw member score writer failed; no PASS");}
					}
				}
				success=true;
			}catch(Throwable caught){failure=caught; next.set(tasks.length);}
		}
		private void flushRaw(){
			assert(rawWriter!=null) : "Only retained-score workers flush raw records";
			if(line.length()==0){return;}
			require(outputBytes.addAndGet(line.length())<=maxOutput,"Raw member score output exceeds byte cap; no PASS");
			rawWriter.print(line);line.clear();
		}
		private void build(Family f){
			assert(f.sequences!=null && f.sequences.length==f.assigned) : "Canonical load must supply every member before profile voting";
			final long start=System.nanoTime();
			final long cpuStart=cpuTime(CPU_BEAN);
			final DDLFamilyProfiles.Profile profile=new DDLFamilyProfiles.Profile(1024);
			final ProteinSpectrum.Spectrum spectrum=new ProteinSpectrum.Spectrum(2);
			for(int i=0;i<f.sequences.length;i++){
				final byte[] sequence=f.sequences[i]; require(sequence!=null,"Missing assigned member before reference construction: "+f.id);
				f.bases+=sequence.length; f.words+=positionMode.positional() ?
					DDLFamilyProfiles.sketchPositional(sequence,geometry,0,partitions,overlap,positionMode,sketch,blacklist) : DDLFamilyProfiles.sketch(sequence,geometry,0,sketch,blacklist);
				profile.add(sketch); spectrum.add(sequence); f.sequences[i]=null;
			}
			require(profile.members()==f.assigned && spectrum.records==f.assigned,"Full reference population must equal assigned count: "+f.id);
			spectrum.validate(); line.clear(); appendProfile(line,f,profile,f.bases,f.words);
			require(outputBytes.addAndGet(line.length())<=maxOutput,"Compact profile output exceeds byte cap; no PASS");
			f.profileRows=line.toBytes(); f.spectrum=spectrum;
			f.buildCpuNanos=cpuElapsed(cpuStart,cpuTime(CPU_BEAN));
			f.buildNanos=System.nanoTime()-start;
			// The voting map belongs only to this invocation, never to the retained Family result.
		}
		private void score(Family f){
			require(f.sequences!=null && f.reference!=null && f.spectrum!=null,"Positive family requires bound members/reference/spectrum: "+f.id);
			final long start=System.nanoTime(),cpuStart=cpuTime(CPU_BEAN); final double[][] scores=new double[3][f.assigned];
			for(int i=0;i<f.sequences.length;i++){
				final byte[] raw=f.sequences[i]; require(raw!=null,"Missing member before positive scoring: "+f.id);
				final ProteinSequence query=buffers.encode(raw);
				querySpectrum.clear(); querySpectrum.add(raw);
				measure(query,querySpectrum,f.reference,f.spectrum,f.norm,metrics);
				scores[0][i]=metrics.cosine; scores[1][i]=metrics.identityPct; scores[2][i]=metrics.blosum;
				if(rawWriter!=null){
					require(f.memberIds!=null && f.memberIds[i]!=null,"Retained score needs its exact source member ID");
					line.append(f.id).tab().append(f.memberIds[i]).tab().append(Double.toString(metrics.identityPct)).tab().append(Double.toString(metrics.blosum)).nl();
					f.memberIds[i]=null;if(line.length()>=65536){flushRaw();}
				}
				f.scored++; f.sequences[i]=null;
			}
			require(f.scored==f.assigned,"All-member scoring population not conserved: "+f.id);
			final long cpuScored=cpuTime(CPU_BEAN); f.scoreCpuNanos=cpuElapsed(cpuStart,cpuScored);
			final long scored=System.nanoTime(); f.scoreNanos=scored-start;
			if(rawWriter!=null){f.oldCutoffs=new Quantile[3];}
			// Preserve the historical CLI choice; preparation callers select retention without a raw writer.
			f.cutoffs=memberCutoffs(scores,rawWriter==null ? MemberRetention.UNIFORM_99 : MemberRetention.METRIC_SPECIFIC,f.oldCutoffs);
			f.quantileCpuNanos=cpuElapsed(cpuScored,cpuTime(CPU_BEAN));
			f.quantileNanos=System.nanoTime()-scored;
		}
		private final Family[] tasks; private final AtomicInteger next; private final AtomicLong outputBytes;
		private final boolean references; private final long maxOutput;
		private final ProteinPositionMode positionMode;
		private final int partitions,overlap;
		private final DDLPentamerBlacklist blacklist;
		private final ProteinKmerGeometry geometry;
		private final Path rawPath;private ByteStreamWriter rawWriter;
		private final char[] sketch=new char[1024]; private final ByteBuilder line=new ByteBuilder();
		private final Buffers buffers=new Buffers(); private final Metrics metrics=new Metrics();
		private final ProteinSpectrum.Spectrum querySpectrum=new ProteinSpectrum.Spectrum(2);
		boolean success; Throwable failure;
	}

	private static final class Family{
		Family(int rank_,String id_,int n){rank=rank_; id=id_; assigned=n;}
		final int rank,assigned; final String id;
		int scored;
		boolean selected; long buildNanos,scoreNanos,quantileNanos,bases,words;
		long buildCpuNanos=-1,scoreCpuNanos=-1,quantileCpuNanos=-1;
		byte[][] sequences,memberIds; byte[] profileRows,consensusRaw; Quantile[] cutoffs,oldCutoffs;
		ProteinSequence reference; ProteinSpectrum.Spectrum spectrum; double norm;
	}

	private static final class Chunk{
		Chunk(Path p,String h,long n,long r){path=p; hash=h; members=n; residues=r;}
		final Path path; final String hash; final long members,residues;
	}

	/** Capacity-only scratch; the owning reader validates its narrower package member range first. */
	static final class Buffers{
		final byte[][] bytes=new byte[201][]; final ProteinSequence[] queries=new ProteinSequence[201];
		byte[] raw(int length){require(length>=40 && length<=200,"Joined-tail scratch supports40..200 residues"); if(bytes[length]==null){bytes[length]=new byte[length];} return bytes[length];}
		void validate(byte[] raw){for(byte b : raw){require(ENCODING[b&255]>=0,"Illegal joined protein residue");}}
		ProteinSequence encode(byte[] raw){
			ProteinSequence q=queries[raw.length];
			if(q==null){q=new ProteinSequence("calibration_query",raw); queries[raw.length]=q;}
			else{for(int i=0;i<raw.length;i++){final int e=ENCODING[raw[i]&255]; require(e>=0,"Illegal joined protein residue"); q.enc[i]=(byte)e;}}
			return q;
		}
	}

	private static byte[] encoding(){
		final byte[] table=new byte[256]; Arrays.fill(table,(byte)-1);
		for(int c=0;c<128;c++){
			try{table[c]=Blosum62.encode(new byte[]{(byte)c},"encoding_table")[0];}catch(RuntimeException unsupported){}
		}
		return table;
	}

	private static Path input(HashMap<String,String> o,String key){final Path p=Paths.get(required(o,key)); pin(p,required(o,key+"sha80")); return p;}
	private static void verifyInput(Path path,String hash,VerificationMode mode,StructuralInputs structural){
		require(mode!=null, "Input verification policy must not be inferred from a missing value");
		if(mode==VerificationMode.SHA80){pin(path,hash);return;}
		require(path!=null && path.isAbsolute() && Files.isRegularFile(path) && !Files.isSymbolicLink(path),
			"Structural verification requires a regular absolute input: "+path);
		require(validIdentity(hash,mode),"Invalid declared structural input identity: "+hash);
		require(structural!=null,"Structural reads require their metadata stability session");
		try{structural.track(path);}
		catch(java.io.IOException | IllegalStateException failure){throw new IllegalArgumentException("Structural input changed while loading: "+path,failure);}
	}
	/** Only explicit structural sessions accept the honest unmeasured identity marker. */
	private static boolean validIdentity(String hash,VerificationMode mode){
		return hash!=null && (hash.matches("[0-9a-f]{20}") || mode==VerificationMode.STRUCTURAL && hash.equals("unmeasured"));
	}
	private static void pin(Path path,String hash){
		require(path.isAbsolute() && Files.isRegularFile(path) && !Files.isSymbolicLink(path),"Require a regular absolute input: "+path);
		require(DigestSuffix.file(path.toString()).equals(DigestSuffix.requireSuffix(hash,"input sha80")),"Input pin mismatch: "+path);
	}
	private static void header(ByteFile input,String value){require(Arrays.equals(input.nextLine(),value.getBytes(StandardCharsets.US_ASCII)),"Unexpected input header: "+value);}
	private static String required(HashMap<String,String> o,String key){final String v=o.get(key); require(v!=null && !v.isEmpty(),"Missing option: "+key); return v;}
	private static long number(HashMap<String,String> o,String key){final long n=Long.parseLong(required(o,key)); require(n>=0,"Negative option: "+key); return n;}
	private static int integer(HashMap<String,String> o,String key){return Math.toIntExact(number(o,key));}
	private static double seconds(long a,long b){return (b-a)*1e-9;}
	private static ByteStreamWriter writer(Path path){final ByteStreamWriter w=new ByteStreamWriter(path.toString(),false,false,false); w.start(); return w;}
	private static void write(Path path,String value) throws Exception{Files.write(path,value.getBytes(StandardCharsets.UTF_8));}
	private static void require(boolean condition,String message){if(!condition){throw new IllegalArgumentException(message);}}
	private static String firstToken(String id){
		require(id!=null && !id.isEmpty(),"Missing shared consensus header");
		int end=0; while(end<id.length() && !Character.isWhitespace(id.charAt(end))){end++;}
		return id.substring(0,end);
	}

	public static final String CUTOFF_HEADER="roster_rank\tfamily_id\tmembers\tstatus\tcosine_cutoff\tcosine_retained\tidentity_pct_cutoff\tidentity_retained\tblosum_a_sq_cutoff\tblosum_retained";
	public static final String RAW_SCORE_HEADER="family_id\tmember_id\tidentity_pct\tblosum_a_sq";
	private static final byte[] ENCODING=encoding();
	private static final ThreadMXBean CPU_BEAN=cpuBean();
	private static final String[] REFERENCE_FILES={"profiles.tsv","spectra.tsv","members.tsv","summary.tsv","recipe.txt","timing.tsv","PASS"};
	private static final java.util.List<String> OPTIONS=Arrays.asList("mode","families","familiessha80","chunks","chunkssha80","models","modelssha80","references","referencessha80",
			"bucket","buckets","familyranks","expectedfamilies","expectedbatches","expectedmembers","maxbucketmembers","maxbucketbytes","maxoutputbytes","maxreferenceaa","t","out","retainscores","consensus","consensussha80",
			"k","aabits","aapartition","sketch","positional","positionmode","partitions","overlap","blacklist","blacklistsha80","minmemberaa","maxmemberaa");
}
