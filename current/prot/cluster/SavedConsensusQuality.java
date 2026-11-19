package prot.cluster;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;

import ddl.DDLPentamerBlacklist;
import ddl.FamilySketchBuilder;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prot.HbmPositionModel;
import prot.ProteinSequence;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.IntHashMap;
import template.Accumulator;
import template.ThreadWaiter;

/** Read-only cohesion and neighbor reports for positive-member families in a saved codec.
 * The complete membership table is validated, including zero-member templates;
 * empty templates remain in the source codec but are absent from both diagnostic rosters.
 */
public final class SavedConsensusQuality{

	public static void main(String[] args) throws Exception{
		final SavedConsensusQuality tool=new SavedConsensusQuality(args);
		try{tool.process();}
		finally{Shared.closeStream(tool.outstream);}
	}

	public SavedConsensusQuality(String[] args){
		final PreParser pre=new PreParser(args,getClass(),false);
		args=pre.args; outstream=pre.outstream;
		final Parser parser=new Parser(); parser.overwrite=true;
		String codec_=null, members_=null, quality_=null;
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			final String key=pair[0].toLowerCase(java.util.Locale.ROOT);
			final String value=pair.length>1 && !pair[1].equalsIgnoreCase("null") ? pair[1] : null;
			if(key.equals("ddl") || key.equals("codec") || key.equals("models")){codec_=value;}
			else if(key.equals("members")){members_=value;}
			else if(key.equals("quality")){quality_=value;}
			else if(key.equals("blacklist")){blacklist=value;}
			else if(key.equals("blacklistsha80")){blacklistSha80=value;}
			else if(key.equals("minspan")){minSpan=Parse.parseIntKMG(value);}
			else if(key.equals("maxalignlen")){maxAlignLength=Parse.parseIntKMG(value);}
			else if(key.equals("pinright")){pinRight=Parse.parseBoolean(value);}
			else if(key.equals("trimleft")){trimLeft=Parse.parseBoolean(value);}
			else if(key.equals("denominator")){
				if("aligned".equalsIgnoreCase(value)){denominator=TailAlignmentGate.Denominator.ALIGNED_QUERY;}
				else if("full".equalsIgnoreCase(value)){denominator=TailAlignmentGate.Denominator.FULL_QUERY;}
				else{throw new IllegalArgumentException("denominator must be aligned or full");}
			}else if(parser.parse(arg,key,value)){}
			else if(Tools.looksLikeInputSequenceStream(arg)){parser.in1=arg;}
			else{throw new IllegalArgumentException("Unknown saved-quality parameter: "+arg);}
		}
		input=parser.in1; codec=codec_; members=members_; quality=quality_;
		overwrite=parser.overwrite; workerThreads=Shared.threads();
		if(input==null || codec==null || members==null || quality==null || workerThreads<1 ||
			parser.in2!=null || parser.out1!=null || parser.out2!=null || parser.outsingle!=null ||
			parser.append || parser.maxReads>=0){
			throw new IllegalArgumentException("Saved quality requires in=, ddl=, members=, quality=, one complete input FASTA and no clustering outputs");
		}
		if(maxAlignLength<1 || minSpan<1 || minSpan>maxAlignLength){
			throw new IllegalArgumentException("Invalid saved-quality alignment gate");
		}
		checkOutputAliases(quality,overwrite,input,codec,members,blacklist);
	}

	private void process() throws Exception{
		Shared.setThreads(1); Shared.AMINO_IN=true; Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=Read.LOWER_CASE_TO_N=Read.DOT_DASH_X_TO_N=Read.FIX_HEADER=false;
		Read.JUNK_MODE=Read.IGNORE_JUNK;
		ReadWrite.setZipThreads(1);
		final ClusteringInput sample=ClusteringInput.load(input,-1,false);
		final ClusterSketchCodec.Loaded loaded=ClusterSketchCodec.read(codec);
		final FamilySketchBuilder sketcher=sketcher(loaded.header);
		final TailAlignmentGate gate=new TailAlignmentGate(pinRight ?
			(trimLeft ? TailAlignmentGate.Mode.PIN_RIGHT_TRIM_LEFT : TailAlignmentGate.Mode.PIN_RIGHT_FULL_QUERY) :
			TailAlignmentGate.Mode.FULL_QUERY_UNPINNED, denominator, maxAlignLength, minSpan, 0, 0);
		write(sample,loaded,members,gate,sketcher,workerThreads,quality,overwrite);
		outstream.println("Measured active-family consensus quality from "+loaded.header.families+
			" saved models and "+sample.sequences.length+" input proteins.");
	}

	public static void write(ClusteringInput sample,ClusterSketchCodec.Loaded loaded,String members,
			TailAlignmentGate gate,FamilySketchBuilder sketcher,int threads,String prefix,
			boolean overwrite){
		if(gate==null || sketcher==null || threads<1){
			throw new IllegalArgumentException("Saved quality needs a gate, sketcher and positive threads");
		}
		checkOutputAliases(prefix,overwrite,members);
		final Family[] families=attachMembers(sample,loaded,members);
		final ClusterSketchCodec.Header header=loaded.header;
		if(!sketcher.matchesGeometry(header.width,header.geometry,header.sketchSeed,header.positionMode,
			header.partitions,header.overlap,header.blacklistSha80)){
			throw new IllegalArgumentException("Saved-quality sketcher differs from the saved codec geometry");
		}
		validateSavedModels(families,gate);
		// Compact only this private attached roster, after every serialized count was checked.
		// An empty template has no empirical cohesion and is not a diagnostic neighbor.
		int count=0;
		for(Family family : families){if(family.size()>0){families[count++]=family;}}
		final Family[] active=count==families.length ? families : Arrays.copyOf(families,count);
		final ConsensusQualityReport cohesion=measure(sample.sequences,active,gate,threads);
		cohesion.write(prefix+".cohesion.tsv",overwrite,active);
		ConsensusNeighborQualityReport.write(active,cohesion,sketcher,gate,threads,prefix,overwrite);
	}

	/** Measures existing nonempty, read-only models without attaching members or rebuilding payloads. */
	public static ConsensusQualityReport measure(ProteinSequence[] sequences,Family[] families,
			TailAlignmentGate gate,int threads){
		if(sequences==null || families==null || gate==null || threads<1){
			throw new IllegalArgumentException("Saved quality needs sequences, families, a gate and positive threads");
		}
		final Work work=new Work(sequences,families,gate,threads);
		work.run();
		final ConsensusQualityReport report=new ConsensusQualityReport();
		for(ConsensusQualityReport.Snapshot snapshot : work.snapshots){report.publish(snapshot);}
		return report;
	}

	private FamilySketchBuilder sketcher(ClusterSketchCodec.Header header) throws Exception{
		if(header.blacklistSha80==null){
			if(blacklist!=null || blacklistSha80!=null){
				throw new IllegalArgumentException("Saved codec declares no blacklist but a blacklist was supplied");
			}
			return new FamilySketchBuilder(header.width,header.geometry,header.sketchSeed,
				header.positionMode,header.partitions,header.overlap,null);
		}
		if(blacklist==null){throw new IllegalArgumentException("Saved codec requires blacklist "+header.blacklistSha80);}
		if(blacklistSha80!=null && !blacklistSha80.equals(header.blacklistSha80)){
			throw new IllegalArgumentException("Configured blacklist sha80 differs from saved codec");
		}
		final DDLPentamerBlacklist loaded=DDLPentamerBlacklist.load(Paths.get(blacklist),header.blacklistSha80);
		return new FamilySketchBuilder(header.width,header.geometry,header.sketchSeed,
			header.positionMode,header.partitions,header.overlap,loaded);
	}

	private static Family[] attachMembers(ClusteringInput sample,ClusterSketchCodec.Loaded loaded,String path){
		if(sample==null || sample.sequences==null || loaded==null || loaded.header==null ||
			loaded.families==null || path==null){
			throw new IllegalArgumentException("Saved quality needs sample, codec and membership inputs");
		}
		if(loaded.header.families!=loaded.families.length){
			throw new IllegalArgumentException("Codec family count disagrees with decoded roster");
		}
		final IntHashMap indexById=new IntHashMap(loaded.families.length*2+1);
		final Family[] families=new Family[loaded.families.length];
		for(int i=0; i<families.length; i++){
			final ClusterSketchCodec.Entry entry=loaded.families[i];
			final Family family=Family.fromSeed(entry,true);
			if(indexById.get(family.familyId)>=0){
				throw new IllegalArgumentException("Duplicate saved family ID: "+family.familyId);
			}
			indexById.put(family.familyId,i);
			families[i]=family;
		}
		final int[] counts=new int[families.length];
		final ByteFile file=ByteFile.makeByteFile1(path,false);
		try{
			final byte[] header=file.nextLine();
			if(header==null || !Arrays.equals(header,MEMBER_HEADER)){
				throw new IllegalArgumentException("Membership TSV must start with the native member header");
			}
			final LineParser1 parser=new LineParser1('\t');
			int expected=0;
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine(), expected++){
				parser.set(line);
				if(parser.terms()!=4 || expected>=sample.sequences.length ||
					parseInt(parser,0,"input_index")!=expected ||
					!parser.termEquals(sample.sequences[expected].id,1)){
					throw new IllegalArgumentException("Membership row does not match input ordinal "+expected);
				}
				final int familyId=parseInt(parser,2,"family_id");
				if(parser.termEquals("assigned",3)){
					final int index=indexById.get(familyId);
					if(familyId<0 || index<0){
						throw new IllegalArgumentException("Assigned row references missing saved family "+familyId);
					}
					families[index].addMember(expected); counts[index]++;
				}else if(parser.termEquals("unassigned",3) || parser.termEquals("filtered_length",3)){
					if(familyId!=-1){throw new IllegalArgumentException("Unassigned membership rows must use family_id=-1");}
				}else{throw new IllegalArgumentException("Unknown membership status at input ordinal "+expected);}
			}
			if(expected!=sample.sequences.length){
				throw new IllegalArgumentException("Membership TSV is incomplete: "+expected+" of "+sample.sequences.length);
			}
		}finally{
			if(file.close()){throw new IllegalStateException("Membership TSV read failed: "+path);}
		}
		for(int i=0; i<families.length; i++){
			final int expected=loaded.families[i].members;
			if(counts[i]!=expected){
				throw new IllegalArgumentException("Membership count for family "+families[i].familyId+
					" differs from saved codec: "+counts[i]+" != "+expected);
			}
		}
		return families;
	}

	private static void validateSavedModels(Family[] families,TailAlignmentGate gate){
		for(Family family : families){
			final byte[] consensus=family.consensusUnsafe();
			if(consensus==null || consensus.length<gate.minSpan){
				throw new IllegalArgumentException("Saved family consensus is shorter than minspan: "+family.familyId);
			}
		}
	}

	private static void checkOutputAliases(String prefix,boolean overwrite,String... protectedSources){
		final String[] neighbors=ConsensusNeighborQualityReport.outputPaths(prefix);
		final ArrayList<String> paths=new ArrayList<String>();
		paths.add(prefix+".cohesion.tsv");
		paths.addAll(Arrays.asList(neighbors));
		final ArrayList<Path> existing=new ArrayList<Path>();
		try{
			for(int i=0; i<paths.size(); i++){
				final String name=paths.get(i);
				if(FileFormat.isStdio(name)){throw new IllegalArgumentException("Quality outputs require disk files: "+name);}
				final Path a=outputPath(name,overwrite);
				for(int j=i+1; j<paths.size(); j++){
					if(samePath(a,outputPath(paths.get(j),overwrite))){
						throw new IllegalArgumentException("Saved-quality outputs alias each other");
					}
				}
				for(String source : protectedSources){
					if(source==null){continue;}
					if(samePath(a,Paths.get(source))){
						throw new IllegalArgumentException("Saved-quality output aliases an input: "+paths.get(i));
					}
				}
				if(Files.exists(a,LinkOption.NOFOLLOW_LINKS)){
					for(Path previous : existing){
						if(Files.isSameFile(previous,a)){throw new IllegalArgumentException("Saved-quality outputs alias each other");}
					}
					existing.add(a);
				}
			}
		}catch(IOException failure){throw new IllegalArgumentException("Could not preflight saved-quality output paths",failure);}
	}

	private static int parseInt(LineParser1 parser,int term,String label){
		final int length=parser.length(term);
		if(length<1){throw new IllegalArgumentException("Empty membership "+label);}
		int i=0, sign=1;
		if(parser.parseByte(term,0)=='-'){
			sign=-1;
			i=1;
			if(i>=length){throw new IllegalArgumentException("Empty membership "+label);}
		}
		long value=0;
		for(; i<length; i++){
			final int digit=parser.parseByte(term,i)-'0';
			if(digit<0 || digit>9){
				throw new IllegalArgumentException("Nondigit membership "+label+": "+field(parser,term));
			}
			value=value*10+digit;
			if(value>Integer.MAX_VALUE+(sign<0 ? 1L : 0L)){
				throw new IllegalArgumentException("Overflowing membership "+label+": "+field(parser,term));
			}
		}
		return sign<0 ? (int)-value : (int)value;
	}

	private static String field(LineParser1 parser,int term){
		final int length=parser.length(term);
		final byte[] bytes=new byte[length];
		for(int i=0; i<length; i++){bytes[i]=parser.parseByte(term,i);}
		return new String(bytes,java.nio.charset.StandardCharsets.US_ASCII);
	}

	private static Path outputPath(String name,boolean overwrite) throws IOException{
		final Path requested=Paths.get(name).toAbsolutePath().normalize();
		final Path parent=requested.getParent();
		if(parent==null || !Files.isDirectory(parent) || !Files.isWritable(parent)){
			throw new IllegalArgumentException("Quality output parent must already be a writable directory: "+name);
		}
		final boolean present=Files.exists(requested,LinkOption.NOFOLLOW_LINKS);
		if(present && (!overwrite || !Files.isRegularFile(requested) || !Files.isWritable(requested))){
			throw new IllegalArgumentException("Quality output exists without overwrite or is not a writable regular file: "+name);
		}
		return present ? requested.toRealPath() : parent.toRealPath().resolve(requested.getFileName());
	}

	private static boolean samePath(Path a,Path b) throws IOException{
		return a.toFile().getCanonicalPath().equals(b.toFile().getCanonicalPath()) ||
			Files.exists(a) && Files.exists(b) && Files.isSameFile(a,b);
	}

	private static final class Work implements Accumulator<Worker>{
		Work(ProteinSequence[] sequences_,Family[] families_,TailAlignmentGate gate_,int threads_){
			sequences=sequences_; families=families_; gate=gate_;
			snapshots=new ConsensusQualityReport.Snapshot[families.length];
			workers=Math.min(threads_,Math.max(1,families.length));
		}
		void run(){
			final ArrayList<Worker> list=new ArrayList<Worker>();
			for(int i=0; i<workers; i++){list.add(new Worker(this));}
			try{
				if(!ThreadWaiter.startAndWait(list,this)){throw new IllegalStateException("Saved quality workers failed");}
			}catch(Throwable cause){
				if(failure==null){failure=cause;}
				final ArrayList<Worker> started=new ArrayList<Worker>();
				for(Worker worker : list){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
				ThreadWaiter.waitForThreadsToFinish(started);
			}
			if(failure!=null){throw new IllegalStateException("Saved quality workers failed",failure);}
			if(done!=families.length){throw new IllegalStateException("Saved quality workers finished "+done+" of "+families.length+" families");}
		}
		int claim(){
			if(failure!=null){return families.length;}
			final int index=next.getAndIncrement();
			return index>=families.length ? families.length : index;
		}
		@Override public void accumulate(Worker worker){
			if(worker.failure!=null && failure==null){failure=worker.failure;}
			done+=worker.done;
		}
		@Override public ReadWriteLock rwlock(){return null;}
		@Override public boolean success(){return failure==null && done==families.length;}
		final ProteinSequence[] sequences;
		final Family[] families;
		final TailAlignmentGate gate;
		final ConsensusQualityReport.Snapshot[] snapshots;
		final AtomicInteger next=new AtomicInteger();
		final int workers;
		int done;
		volatile Throwable failure;
	}

	private static final class Worker extends Thread{
		Worker(Work work_){
			work=work_; scorer=new HbmPositionModel.ConsensusScorer(work.gate.matrix);
			placement=new TailAlignmentGate(work.gate.mode,work.gate.denominator,work.gate.maxLength,1,0,0,null,work.gate.matrix);
		}
		@Override public void run(){
			try{
				for(int i=work.claim(); i<work.families.length; i=work.claim()){
					work.snapshots[i]=score(work.families[i]); done++;
				}
			}catch(Throwable problem){
				failure=problem;
				if(work.failure==null){work.failure=problem;}
			}
		}
		private ConsensusQualityReport.Snapshot score(Family family){
			if(family.size()<1){throw new IllegalStateException("Saved family has no assigned members: "+family.familyId);}
			final byte[] consensus=family.consensusUnsafe();
			scorer.resetConsensus(consensus);
			double identity=0,ratio=0;
			for(int i=0; i<family.size(); i++){
				final int member=family.memberIndex(i);
				if(member<0 || member>=work.sequences.length || work.sequences[member]==null){
					throw new IllegalStateException("Saved family "+family.familyId+" has invalid member "+member);
				}
				final TailAlignmentDecision decision=placement.decide(work.sequences[member].enc,consensus,scorer);
				if(!decision.accepted()){
					if(decision.reason==TailAlignmentDecision.Reason.NO_ALIGNMENT){
						return ConsensusQualityReport.unavailable(family.familyId,family.size(),
							"snapshot_no_positive_placement");
					}
					throw new IllegalStateException("Permissive saved-quality gate rejected member "+member+
						" of family "+family.familyId+": "+decision.reason);
				}
				identity+=decision.identity; ratio+=decision.scoreRatio;
			}
			return ConsensusQualityReport.available(family.familyId,family.size(),identity,ratio);
		}
		final Work work;
		final TailAlignmentGate placement;
		final HbmPositionModel.ConsensusScorer scorer;
		int done;
		Throwable failure;
	}

	private static final byte[] MEMBER_HEADER="input_index\tinput_id\tfamily_id\tstatus".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
	private final PrintStream outstream;
	private final String input,codec,members,quality;
	private String blacklist,blacklistSha80;
	private final boolean overwrite;
	private final int workerThreads;
	private int minSpan=70,maxAlignLength=Integer.MAX_VALUE;
	private boolean pinRight=true,trimLeft=true;
	private TailAlignmentGate.Denominator denominator=TailAlignmentGate.Denominator.ALIGNED_QUERY;
}
