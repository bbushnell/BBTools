package prot;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.concurrent.locks.ReadWriteLock;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.bam.BgzfSettings;
import structures.ByteBuilder;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;

/** Streams accepted winner-member chunks into per-family scalar moments.
 * Adapted from A_SampleByteFileMT: workers parse byte batches, own primitive
 * accumulators, and merge only after joining. Each worker owns an independent
 * chunk reader; an atomic index distributes whole chunks without per-row locks.
 * Shards own disjoint family ranks (rank%shards), including empty roster rows.
 * No train/test filtering, sequence trimming, or biological reassignment occurs.
 * @author Keqing
 */
public final class FamilyMemberScalars implements Accumulator<FamilyMemberScalars.Worker>{

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, FamilyMemberScalars.class, false);
		try{new FamilyMemberScalars(pp.args).process();}finally{Shared.closeStream(pp.outstream);}
	}

	public FamilyMemberScalars(final String[] args) throws Exception{
		final Parser parser=new Parser();
		String minMember=null, maxMember=null;
		for(final String arg : args){
			final String[] kv=Parse.splitOnFirst(arg, '=');
			final String k=kv[0].toLowerCase(Locale.ROOT), v=kv.length>1 ? kv[1] : null;
			if(k.equals("families")){families=v;}
			else if(k.equals("familiessha80")){familyPin=v;}
			else if(k.equals("chunks")){manifest=v;}
			else if(k.equals("chunkssha80")){chunkPin=v;}
			else if(k.equals("expectedfamilies")){nf=Integer.parseInt(v);}
			else if(k.equals("expectedmembers")){expectedMembers=Long.parseLong(v);}
			else if(k.equals("batches")){batches=Integer.parseInt(v);}
			else if(k.equals("buckets")){buckets=Integer.parseInt(v);}
			else if(k.equals("shards")){shards=Integer.parseInt(v);}
			else if(k.equals("shard")){shard=Integer.parseInt(v);}
			else if(k.equals("scalars")){
				require("all".equals(v) || "hydrophobicity".equals(v), "scalars must be all or hydrophobicity");
				hydroOnly="hydrophobicity".equals(v);
			}
			else if(k.equals("verification")){
				require(!verificationSeen && ("sha80".equals(v) || "structural".equals(v)),
					"Specify verification=sha80|structural once");
				verificationSeen=true; structural="structural".equals(v);
			}
			else if(k.equals("minmemberaa")){require(minMember==null && v!=null && !v.isEmpty(), "Duplicate or empty minmemberaa"); minMember=v;}
			else if(k.equals("maxmemberaa")){require(maxMember==null && v!=null && !v.isEmpty(), "Duplicate or empty maxmemberaa"); maxMember=v;}
			else if(!parser.parse(arg, k, v)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		memberRange=MemberLengthRange.from(minMember, maxMember);
		require(families!=null && manifest!=null && parser.out1!=null && !parser.append,
			"Require families=TSV chunks=TSV out=TSV and no append");
		require(nf>0 && expectedMembers>=0 && buckets>0 && batches>0 && shards>0 && buckets%shards==0 &&
			shard>=0 && shard<shards && (long)batches*buckets<=Integer.MAX_VALUE, "Invalid roster/shard dimensions");
		out=parser.out1; overwrite=parser.overwrite;
		final File destination=new File(out).getCanonicalFile();
		require(!destination.equals(new File(families).getCanonicalFile()) && !destination.equals(new File(manifest).getCanonicalFile()),
			"Output must not overwrite an input control");
		require(overwrite || !destination.exists(), "Output exists; use a fresh result path");
		threads=Math.max(1, parser.workers()>0 ? parser.workers() : Shared.threads());
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		// Each compute worker already owns an independent chunk/decoder. A decoder
		// pool per chunk would multiply the requested CPU allocation by worker count.
		BgzfSettings.USE_MULTITHREADED_BGZF=false;
		ids=new String[nf]; expected=new long[nf];
		readRoster(); readManifest(destination);
		final int localFamilies=(nf+shards-1-shard)/shards;
		moments=hydroOnly ? null : new Moments(localFamilies);
		hydro=hydroOnly ? new HydroMeans(localFamilies) : null;
	}

	/** Verifies and validates the entire roster, including unselected ranks. */
	private void readRoster(){
		rosterState=verifyInput(families, familyPin);
		final ByteFile input=ByteFile.makeByteFile1(families, false);
		final LineParser1 lp=new LineParser1('\t');
		final HashSet<String> seen=new HashSet<String>();
		int rank=0; long members=0;
		try{
			require(Arrays.equals(input.nextLine(), ascii("roster_rank\tfamily_id\tassigned")), "Wrong family roster header");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line);
				require(lp.terms()==3 && rank<nf && lp.parseInt(0)==rank, "Family roster must be dense and ordered");
				final String id=lp.parseString(1); final long n=lp.parseLong(2);
				require(!id.isEmpty() && seen.add(id) && n>=0, "Duplicate/empty family identity or negative membership");
				ids[rank]=id; expected[rank]=n; members=Math.addExact(members, n); rank++;
			}
		}finally{require(!input.close(), "Roster read failed");}
		require(rank==nf && members==expectedMembers, "Roster family/member totals differ from expected population");
	}

	/** Retains only this shard's chunks but checks the complete manifest's count algebra. */
	private void readManifest(final File destination) throws Exception{
		manifestState=verifyInput(manifest, chunkPin);
		final ByteFile input=ByteFile.makeByteFile1(manifest, false);
		final LineParser1 lp=new LineParser1('\t');
		final HashSet<String> paths=new HashSet<String>();
		int row=0; long total=0, selected=0;
		try{
			require(Arrays.equals(input.nextLine(), ascii(HbmWinnerPartition.CHUNK_HEADER)), "Wrong chunk manifest header");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line);
				final int bucket=row%buckets, batch=row/buckets;
				require(lp.terms()==6 && row<(long)batches*buckets && lp.parseInt(1)==bucket &&
					lp.termEquals(String.format(Locale.ROOT, "b%02d", batch), 0), "Chunks must be ordered batch then bucket");
				final String path=lp.parseString(2), hash=lp.parseString(3);
				final long n=lp.parseLong(4), residues=lp.parseLong(5);
				require(new File(path).isAbsolute() && paths.add(path) && memberRange.containsTotal(n, residues), "Invalid chunk path/count/member range "+memberRange);
				require(hash.matches("[0-9a-f]{20}") || structural && hash.equals("unmeasured"),
					"Chunk identity must be a lowercase sha80, or unmeasured in structural mode");
				// A fresh output cannot alias any existing chunk. Resolve every path only
				// when overwrite is enabled, avoiding thousands of unnecessary NFS lookups.
				if(overwrite){require(!new File(path).getCanonicalFile().equals(destination), "Output aliases a member chunk");}
				total=Math.addExact(total, n);
				if(bucket%shards==shard){
					chunks.add(new Chunk(path, hash, bucket, n, residues, structural ? new FileState(path) : null));
					selected=Math.addExact(selected, n);
				}
				row++;
			}
		}finally{require(!input.close(), "Chunk manifest read failed");}
		long rosterSelected=0;
		for(int rank=shard; rank<nf; rank+=shards){rosterSelected=Math.addExact(rosterSelected, expected[rank]);}
		require(row==(long)batches*buckets && total==expectedMembers && selected==rosterSelected,
			"Chunk count or global/selected membership conservation failed");
	}

	private void process(){
		final long start=System.nanoTime();
		final ArrayList<Worker> workers=new ArrayList<Worker>();
		for(int i=0; i<Math.min(threads, chunks.size()); i++){workers.add(new Worker(this));}
		final boolean ok=ThreadWaiter.startAndWait(workers, this);
		if(!ok || failure.get()!=null){throw new RuntimeException("Member scalar worker failed; no output is published", failure.get());}
		long members=0;
		for(int rank=shard; rank<nf; rank+=shards){
			final long n=hydroOnly ? hydro.n[rank/shards] : moments.n[rank/shards];
			if(n!=expected[rank]){throw new IllegalStateException("Family membership differs for "+ids[rank]+": "+n+" != "+expected[rank]);}
			members=Math.addExact(members, n);
		}
		verifyStructuralInputs();
		write();
		verifyStructuralInputs();
		System.err.println((hydroOnly ? "FAMILY_MEMBER_HYDROPHOBICITY_PASS" : "FAMILY_MEMBER_SCALARS_PASS")+
			" shard="+shard+" shards="+shards+" families="+(hydroOnly ? hydro.n.length : moments.n.length)+
			" members="+members+" chunks="+chunks.size()+" seconds="+((System.nanoTime()-start)*1e-9));
	}

	@Override public void accumulate(final Worker worker){
		if(!worker.success){failure.compareAndSet(null, new IllegalStateException("Worker did not complete")); return;}
		if(hydroOnly){hydro.merge(worker.hydro);}
		else{moments.merge(worker.stats);}
	}
	@Override public boolean success(){return failure.get()==null;}
	@Override public ReadWriteLock rwlock(){return null;}

	/** Writes one row for every selected roster rank; zero-member statistics are undefined. */
	private void write(){
		if(hydroOnly){writeHydrophobicity(); return;}
		final ByteStreamWriter writer=new ByteStreamWriter(out, overwrite, false, false);
		writer.start();
		try{
			final ByteBuilder bb=new ByteBuilder();
			bb.append("family_rank\tfamily_id\tN");
			for(final String name : NAMES){for(final String metric : METRICS){bb.tab().append(name).append('_').append(metric);}}
			writer.print(bb.nl());
			for(int rank=shard; rank<nf; rank+=shards){
				final int family=rank/shards; final long n=moments.n[family];
				bb.clear().append(rank).tab().append(ids[rank]).tab().append(n);
				for(int scalar=0; scalar<NAMES.length; scalar++){
					final int i=family*NAMES.length+scalar;
					if(n==0){bb.append("\tNA\tNA\tNA\tNA");}
					else{
						require(moments.m2[i]>=0 && Double.isFinite(moments.m2[i]), "Nonfinite or negative scalar variance");
						bb.tab().append(moments.mean[i], 12).tab().append(Math.sqrt(moments.m2[i]/n), 12)
							.tab().append(moments.min[i], 12).tab().append(moments.max[i], 12);
					}
				}
				writer.print(bb.nl());
			}
		}finally{require(!writer.poisonAndWait(), "Member scalar output failed");}
	}

	/** Equal weight per retained member, not a pooled-residue mean or consensus substitute. */
	private void writeHydrophobicity(){
		final ByteStreamWriter writer=new ByteStreamWriter(out, overwrite, false, false); writer.start();
		try{
			writer.print("family_rank\tfamily_id\tN\thydrophobicity_mean\n");
			final ByteBuilder bb=new ByteBuilder();
			for(int rank=shard; rank<nf; rank+=shards){
				final int at=rank/shards;final long n=hydro.n[at];
				bb.clear().append(rank).tab().append(ids[rank]).tab().append(n).tab();
				if(n==0){bb.append("NA");}
				else{
					require(Double.isFinite(hydro.mean[at]) && hydro.mean[at]>=-4.5 && hydro.mean[at]<=4.5,
						"Member KD mean must stay within the native residue scale for "+ids[rank]);
					bb.append(hydro.mean[at],12);
				}
				writer.print(bb.nl());
			}
		}finally{require(!writer.poisonAndWait(), "Member hydrophobicity output failed");}
	}

	/** Computes only KD from the existing native table; no pI or other scalar work. */
	private static double memberHydrophobicity(final byte[] sequence,final int length){
		require(length>0 && length<=sequence.length,"Member KD requires a complete nonempty peptide");
		double sum=0;
		for(int i=0;i<length;i++){
			final int aa=ConsensusScalars.aaIndex(sequence[i]);
			require(aa>=0,"Accepted original tail contains a noncanonical amino acid");
			sum+=ConsensusScalars.KD[aa];
		}
		return sum/length;
	}

	/** Independent readers mean worker exceptions cannot strand a shared batch producer. */
	static final class Worker extends Thread{
		Worker(final FamilyMemberScalars owner_){
			owner=owner_;
			stats=owner.hydroOnly ? null : new Moments(owner.moments.n.length);
			hydro=owner.hydroOnly ? new HydroMeans(owner.hydro.n.length) : null;
		}
		@Override public void run(){
			try{
				for(int index=owner.next.getAndIncrement(); index<owner.chunks.size() && owner.failure.get()==null;
						index=owner.next.getAndIncrement()){read(owner.chunks.get(index));}
				success=owner.failure.get()==null;
			}catch(Throwable error){owner.failure.compareAndSet(null, error);}
		}
		private void read(final Chunk chunk){
			if(chunk.state==null){pin(chunk.path, chunk.hash);}else{chunk.state.verify();}
			final File file=new File(chunk.path); final long size=file.length(), modified=file.lastModified();
			final ByteFile input=ByteFile.makeByteFile1(chunk.path, false);
			long members=0, residues=0;
			try{
				require(Arrays.equals(input.nextLine(), ascii(HbmWinnerPartition.MEMBER_HEADER)), "Wrong joined-member header");
				for(ListNum<byte[]> list=input.nextList(); list!=null; list=input.nextList()){
					if(owner.failure.get()!=null){return;}
					for(final byte[] line : list){
						lp.set(line); require(lp.terms()==6, "Joined member must have six columns");
						final int rank=lp.parseInt(0), len=lp.length(5);
						require(rank>=0 && rank<owner.nf && rank%owner.buckets==chunk.bucket &&
							lp.termEquals(owner.ids[rank], 1), "Member family rank/identity/bucket mismatch");
						require(owner.memberRange.contains(len) && lp.length(2)>0 && lp.length(3)>0 && lp.parseLong(4)>=0,
							"Invalid original tail length/identity/ordinal");
						lp.setBounds(5); System.arraycopy(line, lp.a(), sequence, 0, len);
						if(owner.hydroOnly){hydro.add(rank/owner.shards,memberHydrophobicity(sequence,len));}
						else{
							final int valid=ConsensusScalars.computeScalars(sequence, len, scalars, 0,scalarMask);
							require(valid==len, "Accepted original tail contains a noncanonical amino acid");
							require(scalarMask[0]==ConsensusScalars.ALL_DEFINED,"Member statistics require all five measured scalars");
							stats.add(rank/owner.shards, scalars);
						}
						members++; residues+=len;
					}
				}
			}finally{require(!input.close(), "Member chunk read failed");}
			if(members!=chunk.members || residues!=chunk.residues){throw new IllegalStateException("Chunk counts differ: "+chunk.path);}
			require(file.length()==size && file.lastModified()==modified, "Pinned chunk changed while reading");
			if(chunk.state!=null){chunk.state.verify();}
		}
		final FamilyMemberScalars owner;
		final Moments stats;
		final HydroMeans hydro;
		final LineParser1 lp=new LineParser1('\t');
		final byte[] sequence=new byte[200];
		final double[] scalars=new double[5];
		final int[] scalarMask=new int[1];
		boolean success;
	}

	/** Primitive online member means and coarse worker reduction; no per-record objects or unused scalar fields. */
	private static final class HydroMeans{
		HydroMeans(final int families){
			require(families>=0,"Hydrophobicity accumulator needs a valid local roster size");
			n=new long[families];mean=new double[families];
		}
		void add(final int family,final double value){
			require(family>=0 && family<n.length && Double.isFinite(value) && value>=-4.5 && value<=4.5,
				"Every retained member contributes one finite KD mean to its own family");
			final long count=Math.addExact(n[family],1);n[family]=count;
			mean[family]+=(value-mean[family])/count;
		}
		void merge(final HydroMeans other){
			require(other.n.length==n.length,"Cannot combine different hydrophobicity rosters");
			for(int i=0;i<n.length;i++){
				final long b=other.n[i];if(b==0){continue;}
				final long total=Math.addExact(n[i],b);
				mean[i]+=(other.mean[i]-mean[i])*((double)b/total);n[i]=total;
			}
		}
		final long[] n;
		final double[] mean;
	}

	/** Welford updates and parallel moment combination, with no per-member allocation. */
	public static final class Moments{
		public Moments(final int families){
			require(families>=0 && families<=Integer.MAX_VALUE/5, "Invalid moments dimensions");
			n=new long[families]; mean=new double[families*5]; m2=new double[mean.length];
			min=new double[mean.length]; max=new double[mean.length];
			Arrays.fill(min, Double.POSITIVE_INFINITY); Arrays.fill(max, Double.NEGATIVE_INFINITY);
		}
		public void add(final int family, final double[] values){
			assert(family>=0 && family<n.length && values.length>=5) : "Scalar accumulator needs a valid local family and five values";
			final long count=Math.addExact(n[family], 1); n[family]=count;
			for(int s=0; s<5; s++){
				final double value=values[s]; require(Double.isFinite(value), "Member scalar is not finite");
				final int i=family*5+s; final double delta=value-mean[i];
				mean[i]+=delta/count; m2[i]+=delta*(value-mean[i]);
				min[i]=Math.min(min[i], value); max[i]=Math.max(max[i], value);
			}
		}
		public void merge(final Moments other){
			require(other.n.length==n.length, "Cannot merge different family dimensions");
			for(int f=0; f<n.length; f++){
				final long a=n[f], b=other.n[f]; if(b==0){continue;}
				final long total=Math.addExact(a, b);
				for(int s=0; s<5; s++){
					final int i=f*5+s; final double delta=other.mean[i]-mean[i];
					m2[i]+=other.m2[i]+delta*delta*((double)a*b/total);
					mean[i]+=delta*((double)b/total);
					min[i]=Math.min(min[i], other.min[i]); max[i]=Math.max(max[i], other.max[i]);
				}
				n[f]=total;
			}
		}
		/** Exports five mean/SD/min/max groups through the existing 12-decimal table boundary.
		 * Input moments are unchanged; out has exactly twenty slots and scratch is cleared/reused.
		 * Empty families have no measured statistics and must be handled separately by the caller.
		 * This deliberately uses the same formatter and correctly rounded parser as file export.
		 */
		public void exportStatistics(final int family,final double[] out,final ByteBuilder scratch){
			require(family>=0 && family<n.length && n[family]>0,"Scalar statistics require a nonempty family in this accumulator");
			require(out!=null && out.length==20 && scratch!=null,"Scalar export needs twenty output slots and caller-owned formatting scratch");
			require(out!=mean && out!=m2 && out!=min && out!=max,"Scalar output must not alias accumulated moments");
			for(int s=0;s<5;s++){
				final int i=family*5+s,offset=s*4;
				require(m2[i]>=0 && Double.isFinite(m2[i]),"Nonfinite or negative scalar variance");
				out[offset]=roundedStatistic(mean[i],scratch);
				out[offset+1]=roundedStatistic(Math.sqrt(m2[i]/n[family]),scratch);
				out[offset+2]=roundedStatistic(min[i],scratch);
				out[offset+3]=roundedStatistic(max[i],scratch);
			}
		}
		private static double roundedStatistic(final double value,final ByteBuilder scratch){
			require(Double.isFinite(value),"Only measured finite scalar statistics can cross the table boundary");
			return Double.parseDouble(scratch.clear().append(value,12).toString());
		}
		public final long[] n;
		public final double[] mean, m2, min, max;
	}

	private static final class Chunk{
		Chunk(String path_, String hash_, int bucket_, long members_, long residues_, FileState state_){
			path=path_; hash=hash_; bucket=bucket_; members=members_; residues=residues_; state=state_;
		}
		final String path, hash; final int bucket; final long members, residues;
		final FileState state;
	}
	/** Exact mode retains its original checks; structural mode binds file metadata instead. */
	private FileState verifyInput(final String path, final String identity){
		if(!structural){pin(path, identity);return null;}
		require(identity!=null && (identity.equals("unmeasured") || identity.matches("[0-9a-f]{20}")),
			"Invalid structural member input identity: "+path);
		return new FileState(path);
	}
	/** The whole declared input view must remain stable until the final result is written. */
	private void verifyStructuralInputs(){
		if(!structural){return;}
		rosterState.verify();manifestState.verify();
		for(Chunk chunk:chunks){chunk.state.verify();}
	}
	/** Immutable per-file metadata captured before any worker can read it; never a content hash. */
	private static final class FileState{
		FileState(final String name){path=Paths.get(name);initial=read();}
		private BasicFileAttributes read(){
			require(path.isAbsolute() && !Files.isSymbolicLink(path), "Structural member input must be absolute and nonsymlink: "+path);
			try{
				final BasicFileAttributes state=Files.readAttributes(path, BasicFileAttributes.class, LinkOption.NOFOLLOW_LINKS);
				require(state.isRegularFile(), "Structural member input must be a regular file: "+path);
				return state;
			}catch(java.io.IOException e){throw new IllegalArgumentException("Cannot inspect structural member input: "+path,e);}
		}
		void verify(){
			final BasicFileAttributes current=read();
			require(initial.size()==current.size() && initial.lastModifiedTime().equals(current.lastModifiedTime()) &&
				java.util.Objects.equals(initial.fileKey(),current.fileKey()), "Structural member input metadata changed: "+path);
		}
		final Path path;
		final BasicFileAttributes initial;
	}
	private static void pin(final String path, final String hash){
		require(hash!=null && hash.matches("[0-9a-f]{20}"), "Input digest must be a lowercase sha80");
		require(DigestSuffix.file(path).equals(hash), "Input content differs: "+path);
	}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}
	private static byte[] ascii(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static final String[] NAMES={"hydrophobicity","net_charge","aromaticity","fraction_small","pI_tail"};
	private static final String[] METRICS={"mean","sd","min","max"};
	private String families, manifest, out, familyPin, chunkPin;
	private final MemberLengthRange memberRange;
	private int nf=79565, buckets=256, batches=64, shards=256, shard=-1, threads;
	private long expectedMembers=124056489;
	private boolean overwrite,hydroOnly,structural=true,verificationSeen;
	private FileState rosterState,manifestState;
	private String[] ids;
	private long[] expected;
	private Moments moments;
	private HydroMeans hydro;
	private final ArrayList<Chunk> chunks=new ArrayList<Chunk>();
	private final AtomicInteger next=new AtomicInteger();
	private final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
}
