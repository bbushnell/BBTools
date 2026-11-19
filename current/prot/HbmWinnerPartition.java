package prot;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.zip.GZIPInputStream;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.bam.BgzfOutputStream;
import structures.ByteBuilder;
import structures.ListNum;

/** Joins frozen winner assignments to their original sequences, without rescoring.
 * One process handles one batch in manifest order; independent batches may run in
 * parallel. Scalar BGZF writers avoid a compression thread pool per bucket.
 * @author Keqing
 */
public final class HbmWinnerPartition {

	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args, HbmWinnerPartition.class, false);
			try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	/** Throwing production entry point, also used by finite file-backed fixtures. */
	public static void run(String[] args) throws Exception{
		final HashMap<String,String> o=new HashMap<String,String>();
		for(String arg : args){
			final String[] kv=Parse.splitOnFirst(arg, '=');
			require(kv.length==2, "Expected flag=value");
			final String key=kv[0].toLowerCase(Locale.ROOT);
			require(Arrays.asList("mode","in","insha80","models","modelssha80","support","supportsha80",
				"expectedqueries","expectedassigned","expectedgenomes","buckets","batch","parts","out","maxbytes").contains(key), "Unknown option: "+key);
			require(o.put(key, kv[1])==null, "Duplicate option: "+key);
		}
		Shared.setThreads(1); Shared.AMINO_IN=true;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final Inputs inputs=new Inputs(o);
		final Path out=Paths.get(required(o, "out")).toAbsolutePath().normalize();
		require(!Files.exists(out, LinkOption.NOFOLLOW_LINKS) && Files.isDirectory(out.getParent()), "Output must be fresh under an existing parent");
		final String mode=required(o, "mode");
		if(mode.equals("partition")){partition(inputs, required(o, "batch"), out, number(required(o, "maxbytes")));}
		else if(mode.equals("collect")){collect(inputs, Paths.get(required(o, "parts")).toAbsolutePath().normalize(), out);}
		else{throw new IllegalArgumentException("mode must be partition or collect");}
	}

	private static void partition(Inputs inputs,String batch,Path out,long maxBytes) throws Exception{
		require(maxBytes>0, "Positive decoded-output cap required");
		final String[] binding=inputs.batch(batch);
		final Batch b=new Batch(inputs, binding);
		Files.createDirectory(out);
		final Chunk[] chunks=new Chunk[inputs.buckets];
		final long[] counts=new long[inputs.ids.length], total=new long[8];
		long bytes=0;
		try{
			for(int i=0; i<chunks.length; i++){
				chunks[i]=new Chunk(out.resolve(String.format(Locale.ROOT, "bucket_%03d.tsv.gz", i)));
				bytes=Math.addExact(bytes, chunks[i].bytes);
			}
			require(bytes<=maxBytes, "Bucket headers exceed decoded-output cap");
			final HashSet<String> genomes=new HashSet<String>();
			for(String[] shard : b.shards){
				require(safe(shard[0]) && genomes.add(shard[0]), "Unsafe or repeated genome in batch");
				final long[] got=join(inputs, b, shard, chunks, counts, maxBytes-bytes);
				for(int i=0; i<8; i++){total[i]=Math.addExact(total[i], got[i]);}
				bytes=Math.addExact(bytes, got[8]);
				System.err.println("winner_partition_genome\t"+batch+"\t"+shard[0]+"\t"+total[0]+"\t"+total[2]);
			}
			require(Arrays.equals(total, b.summary), "Batch source/assignment/summary totals differ");
		}finally{
			Exception failure=null;
			for(Chunk chunk : chunks){
				if(chunk!=null){try{chunk.close();}catch(Exception e){if(failure==null){failure=e;}else{failure.addSuppressed(e);}}}
			}
			if(failure!=null){throw failure;}
		}
		final StringBuilder manifest=new StringBuilder(CHUNK_HEADER).append('\n');
		final StringBuilder verification=new StringBuilder("bucket\tdecoded_bytes\tdecoded_sha80\n");
		for(int i=0; i<chunks.length; i++){
			final Chunk chunk=chunks[i]; chunk.verify();
			manifest.append(batch).append('\t').append(i).append('\t').append(chunk.path).append('\t')
				.append(DigestSuffix.file(chunk.path.toString())).append('\t').append(chunk.members).append('\t').append(chunk.residues).append('\n');
			verification.append(i).append('\t').append(chunk.bytes).append('\t').append(chunk.decodedPin).append('\n');
		}
		write(out.resolve("chunks.tsv"), manifest.toString());
		write(out.resolve("compression.tsv"), verification.toString());
		writeCounts(out.resolve("families.tsv"), inputs.ids, counts);
		write(out.resolve("summary.tsv"), SUMMARY_HEADER+"\n"+numbers(total)+"\n");
		write(out.resolve("recipe.tsv"), inputs.recipe()+"batch\t"+batch+"\nbatch_recipe_sha80\t"+binding[4]+"\nmax_decoded_bytes\t"+maxBytes+"\n");
		b.verify(); inputs.verify();
		publish(out, "HBM_WINNER_PARTITION_PASS", new String[]{"chunks.tsv","families.tsv","summary.tsv","recipe.tsv","compression.tsv"});
		System.err.println("HBM_WINNER_PARTITION_PASS batch="+batch+" queries="+total[0]+" assigned="+total[2]);
	}

	/** Source and calls advance in lockstep, including rejected records; no trimming or retranslation. */
	private static long[] join(Inputs inputs,Batch batch,String[] s,Chunk[] chunks,long[] family,long byteCap) throws Exception{
		assert(chunks.length==inputs.buckets && family.length==inputs.ids.length) : "Bucket and family arrays use the pinned roster";
		final String genome=s[0], prefix=genome+"|";
		final long nq=number(s[1]), nr=number(s[2]);
		final Path recorded=Paths.get(s[3]);
		require(recorded.isAbsolute() && recorded.normalize().equals(recorded),"Original FASTA path must be absolute and normalized");
		final Path fasta=regular(recorded.toRealPath()), dir=batch.root.resolve("genomes").resolve(genome);
		pin(fasta, s[4]);
		final Path complete=dir.resolve("complete.tsv");
		final String completionPin=DigestSuffix.file(regular(complete).toString());
		final HashMap<String,String> c=kv(complete, completionPin);
		require(c.size()==8 && required(c,"schema").equals("winner_shard_v1") && required(c,"recipe_sha80").equals(batch.binding[4])
			&& required(c,"genome").equals(genome) && required(c,"input_sha80").equals(s[4])
			&& number(required(c,"queries"))==nq && number(required(c,"residues"))==nr, "Completion receipt does not bind this source genome");
		final Path calls=dir.resolve("assignments.tsv"), summary=dir.resolve("summary.tsv");
		final long[] saved=summary(summary, required(c,"summary_sha80"));
		require(saved[0]==nq && saved[1]==nr, "Genome manifest/summary mismatch");
		pin(calls, required(c,"assignments_sha80"));
		final ByteFile callReader=ByteFile.makeByteFile1(calls.toString(), false);
		final Streamer source=StreamerFactory.makeStreamer(FileFormat.testInput(fasta.toString(),FileFormat.FASTA,null,true,false),null,true,-1,false,true,0);
		final LineParser1 lp=new LineParser1('\t');
		final HashSet<String> seen=new HashSet<String>();
		final long[] observed=new long[9]; final ByteBuilder row=new ByteBuilder();
		try{
			require(Arrays.equals(callReader.nextLine(), HbmWinnerCalibration.ASSIGNMENT_HEADER.getBytes(StandardCharsets.US_ASCII)), "Expected21-column winner assignments");
			source.start();
			for(ListNum<Read> list=source.nextList(); list!=null; list=source.nextList()){
				for(Read read : list){
					require(read!=null && read.mate==null && read.bases!=null && observed[0]<nq, "Extra or invalid source record");
					final String id=HbmCompetitiveAssign.queryId(read.id);
					require(token(id) && id.startsWith(prefix) && id.length()>prefix.length() && seen.add(id), "Duplicate or wrong-namespace source ID");
					require(read.length()>=100 && read.length()<=200, "Winner input must be an unmodified100..200AA tail");
					for(byte residue : read.bases){
						final int aa=residue>=0 && residue<AminoAcid.acidToNumberNoStops.length ? AminoAcid.acidToNumberNoStops[residue] : -1;
						require(aa>=0 && aa<20, "Ambiguous or stop residue in frozen winner input");
					}
					final byte[] call=callReader.nextLine(); require(call!=null, "Assignment stream ended before source"); lp.set(call);
					require(lp.terms()==21 && lp.termEquals(id,0), "Assignment/source ID or field count differs");
					final int rank=integer(lp,2), mask=integer(lp,20), ties=integer(lp,5), scored=integer(lp,6), nonpositive=integer(lp,7);
					require(rank>=0 && rank<inputs.ids.length && lp.termEquals(inputs.ids[rank],3) && mask>=0 && mask<=15
						&& lp.termEquals(mask==0 ? "ASSIGNED" : "REJECTED",1) && ties>=1 && ties<=200 && scored==200 && nonpositive>=0 && nonpositive<=200,
						"Winner rank/identity/status/work counters differ");
					if(mask==0){
						row.clear().append(rank).tab().append(inputs.ids[rank]).tab().append(id).tab().append(genome).tab().append(observed[0]).tab().append(read.bases).nl();
						observed[8]=Math.addExact(observed[8],row.length()); require(observed[8]<=byteCap, "Decoded partition output cap exceeded");
						chunks[rank%chunks.length].add(row,read.length()); family[rank]=Math.addExact(family[rank],1); observed[2]++;
					}else{observed[3]++;}
					observed[0]++; observed[1]=Math.addExact(observed[1],read.length()); observed[4]+=scored; observed[5]+=nonpositive;
					if(ties>1){observed[6]++;} observed[7]++;
				}
			}
			require(callReader.nextLine()==null, "Extra assignments after source EOF");
		}finally{
			try{source.close(); require(!source.errorState(), "Source reader failed");}
			finally{require(!callReader.close(), "Assignment reader failed");}
		}
		require(Arrays.equals(Arrays.copyOf(observed,8),saved), "Joined genome count/residue conservation failed");
		pin(fasta,s[4]); require(recorded.toRealPath().equals(fasta),"Original FASTA alias changed during join");
		pin(calls,required(c,"assignments_sha80")); pin(summary,required(c,"summary_sha80")); pin(complete,completionPin);
		return observed;
	}

	/** Metadata-only aggregation; consumers recheck stored chunk pins before parsing. */
	private static void collect(Inputs inputs,Path parts,Path out) throws Exception{
		final long[] counts=new long[inputs.ids.length], total=new long[8];
		final StringBuilder manifest=new StringBuilder(CHUNK_HEADER).append('\n');
		final HashSet<String> genomes=new HashSet<String>();
		for(String[] binding : inputs.batches){
			final Batch batch=new Batch(inputs,binding); final Path dir=parts.resolve(binding[0]);
			for(String[] shard : batch.shards){require(genomes.add(shard[0]), "Genome repeated across batches");}
			final HashMap<String,String> artifacts=readPass(dir,"HBM_WINNER_PARTITION_PASS");
			final HashMap<String,String> recipe=kv(dir.resolve("recipe.tsv"),required(artifacts,"recipe.tsv"));
			for(String key : new String[]{"in_sha80","models_sha80","support_sha80","buckets"}){
				require(required(recipe,key).equals(inputs.recipeValues().get(key)), "Mixed partition input or bucket recipe: "+key);
			}
			require(required(recipe,"batch").equals(binding[0]) && required(recipe,"batch_recipe_sha80").equals(binding[4]), "Partition belongs to another batch");
			final ArrayList<String[]> familyRows=table(dir.resolve("families.tsv"),required(artifacts,"families.tsv"),FAMILY_HEADER,3);
			require(familyRows.size()==counts.length, "Incomplete partition family inventory");
			long assigned=0;
			for(int i=0; i<counts.length; i++){
				final String[] r=familyRows.get(i); require(number(r[0])==i && r[1].equals(inputs.ids[i]), "Partition family order/identity mismatch");
				final long n=number(r[2]); counts[i]=Math.addExact(counts[i],n); assigned=Math.addExact(assigned,n);
			}
			final long[] local=summary(dir.resolve("summary.tsv"),required(artifacts,"summary.tsv"));
			require(Arrays.equals(local,batch.summary) && assigned==local[2], "Partition summary differs from frozen batch");
			for(int i=0; i<8; i++){total[i]=Math.addExact(total[i],local[i]);}
			final ArrayList<String[]> chunkRows=table(dir.resolve("chunks.tsv"),required(artifacts,"chunks.tsv"),CHUNK_HEADER,6);
			require(chunkRows.size()==inputs.buckets, "Missing bucket metadata");
			long chunkTotal=0;
			for(int i=0; i<inputs.buckets; i++){
				final String[] r=chunkRows.get(i);
				require(r[0].equals(binding[0]) && number(r[1])==i && Paths.get(r[2]).equals(dir.resolve(String.format(Locale.ROOT,"bucket_%03d.tsv.gz",i))), "Wrong bucket path/order");
				regular(Paths.get(r[2])); DigestSuffix.requireSuffix(r[3],"chunk pin");
				final long n=number(r[4]), residues=number(r[5]);
				require(residues>=Math.multiplyExact(n,100) && residues<=Math.multiplyExact(n,200), "Invalid bucket residues");
				chunkTotal=Math.addExact(chunkTotal,n); manifest.append(String.join("\t",r)).append('\n');
			}
			require(chunkTotal==assigned, "Bucket counts do not conserve batch assignments"); batch.verify();
		}
		require(genomes.size()==inputs.expectedGenomes && total[0]==inputs.expectedQueries && total[2]==inputs.expectedAssigned
			&& Arrays.equals(counts,inputs.assigned), "Global assigned-member census differs from accepted oracle");
		inputs.verify(); Files.createDirectory(out);
		write(out.resolve("chunks.tsv"),manifest.toString()); writeCounts(out.resolve("families.tsv"),inputs.ids,counts);
		write(out.resolve("summary.tsv"),SUMMARY_HEADER+"\n"+numbers(total)+"\n");
		write(out.resolve("recipe.tsv"),inputs.recipe()+"chunk_verification\tproducer_roundtrip_then_consumer_stored_pin\n");
		publish(out,"HBM_WINNER_PARTITION_COLLECTION_PASS",new String[]{"chunks.tsv","families.tsv","summary.tsv","recipe.tsv"});
		System.err.println("HBM_WINNER_PARTITION_COLLECTION_PASS assigned="+total[2]);
	}

	private static final class Inputs {
		Inputs(HashMap<String,String> options) throws Exception{
			o=options; buckets=Math.toIntExact(number(required(o,"buckets"))); require(buckets>0 && buckets<=256,"Buckets must be1..256");
			expectedQueries=number(required(o,"expectedqueries")); expectedAssigned=number(required(o,"expectedassigned"));
			expectedGenomes=number(required(o,"expectedgenomes")); require(expectedAssigned<=expectedQueries && expectedGenomes>0,"Invalid population bounds");
			batches=table(Paths.get(required(o,"in")),required(o,"insha80"),BATCH_HEADER,6);
			require(!batches.isEmpty(),"Empty batch manifest");
			for(int i=0; i<batches.size(); i++){require(batches.get(i)[0].equals(String.format(Locale.ROOT,"b%02d",i)),"Batch IDs must be dense and ordered");}
			final ArrayList<String[]> models=table(Paths.get(required(o,"models")),required(o,"modelssha80"),HbmTailLibrary.MODELS_HEADER,9);
			final ArrayList<String[]> support=table(Paths.get(required(o,"support")),required(o,"supportsha80"),SUPPORT_HEADER,5);
			require(!models.isEmpty() && support.size()==models.size(),"Model/support roster size differs");
			ids=new String[models.size()]; assigned=new long[ids.length]; final HashSet<String> seen=new HashSet<String>(); long sum=0;
			for(int i=0; i<ids.length; i++){
				final String[] m=models.get(i), s=support.get(i);
				require(number(m[0])==i && number(s[0])==i && safe(m[1]) && seen.add(m[1]) && m[1].equals(s[1]) && number(m[3])==number(s[2]),"Model/support rank or identity mismatch");
				ids[i]=m[1]; assigned[i]=number(s[3]); sum=Math.addExact(sum,assigned[i]);
				require(Long.parseLong(s[4])==assigned[i]-number(m[3]),"Support delta does not conserve membership");
			}
			require(sum==expectedAssigned,"Support oracle differs from expected assigned population");
		}
		String[] batch(String id){for(String[] b : batches){if(b[0].equals(id)){return b;}} throw new IllegalArgumentException("Unknown batch: "+id);}
		void verify() throws Exception{for(String key : new String[]{"in","models","support"}){pin(Paths.get(required(o,key)),required(o,key+"sha80"));}}
		HashMap<String,String> recipeValues(){
			final HashMap<String,String> r=new HashMap<String,String>();
			for(String key : new String[]{"in","models","support"}){r.put(key+"_sha80",required(o,key+"sha80"));}
			r.put("buckets",Integer.toString(buckets)); return r;
		}
		String recipe(){return "schema\twinner_member_partition_v1\nin_sha80\t"+required(o,"insha80")+"\nmodels_sha80\t"+required(o,"modelssha80")
			+"\nsupport_sha80\t"+required(o,"supportsha80")+"\nbuckets\t"+buckets+"\norder\tbatch_manifest_genome_manifest_source_ordinal\ncompression\tBGZF9_scalar\n";}
		final HashMap<String,String> o; final int buckets; final long expectedQueries,expectedAssigned,expectedGenomes;
		final ArrayList<String[]> batches; final String[] ids; final long[] assigned;
	}

	private static final class Batch {
		Batch(Inputs inputs,String[] row) throws Exception{
			binding=row; root=Paths.get(row[3]); require(root.isAbsolute() && root.normalize().equals(root) && root.toRealPath().equals(root),"Noncanonical result root");
			final HashMap<String,String> r=kv(root.resolve("recipe.tsv"),row[4]);
			require(required(r,"schema").equals("winner_batch_v1") && required(r,"models_sha80").equals(required(inputs.o,"modelssha80"))
				&& required(r,"models").equals(required(inputs.o,"models")) && required(r,"shards").equals(row[1]) && required(r,"shards_sha80").equals(row[2])
				&& required(r,"patience").equals("0") && required(r,"shortlist").equals("200") && required(r,"hbm_kernel").equals("alignTail")
				&& required(r,"policy").equals("hbm_winner_then_blosum_no_fallback"),"Wrong winner policy or recipe binding");
			shards=table(Paths.get(row[1]),row[2],SHARD_HEADER,5); summary=summary(root.resolve("summary.tsv"),row[5]);
			long q=0,n=0; for(String[] s : shards){q=Math.addExact(q,number(s[1])); n=Math.addExact(n,number(s[2]));}
			require(q==summary[0] && n==summary[1],"Batch manifest population mismatch"); verify();
		}
		void verify() throws Exception{
			pin(root.resolve("recipe.tsv"),binding[4]); pin(root.resolve("summary.tsv"),binding[5]); pin(Paths.get(binding[1]),binding[2]);
			require(!Files.exists(root.resolve("RUNNING"),LinkOption.NOFOLLOW_LINKS),"Winner batch is still running");
			require(new String(Files.readAllBytes(regular(root.resolve("PASS"))),StandardCharsets.US_ASCII).equals(
				"HBM_WINNER_BATCH_PASS\nrecipe_sha80\t"+binding[4]+"\nsummary_sha80\t"+binding[5]+"\n"),"Winner batch lacks matching PASS");
		}
		final String[] binding; final Path root; final ArrayList<String[]> shards; final long[] summary;
	}

	private static final class Chunk {
		Chunk(Path path_) throws Exception{
			path=path_; digest=MessageDigest.getInstance("SHA-256");
			out=new BgzfOutputStream(new BufferedOutputStream(Files.newOutputStream(path),65536),9);
			final byte[] header=(MEMBER_HEADER+"\n").getBytes(StandardCharsets.US_ASCII);
			out.write(header); digest.update(header); bytes=header.length;
		}
		void add(ByteBuilder row,int length) throws Exception{
			assert(length>=100 && length<=200) : "Validated source tails alone enter the assigned-member partition";
			out.write(row.array,0,row.length()); digest.update(row.array,0,row.length());
			bytes=Math.addExact(bytes,row.length()); members++; residues+=length;
		}
		void close() throws Exception{out.flush(); out.writeEOF(); out.close(); decodedPin=DigestSuffix.fromDigest(digest.digest());}
		void verify() throws Exception{
			final MessageDigest decoded=MessageDigest.getInstance("SHA-256"); final byte[] buffer=new byte[65536]; long n=0;
			try(InputStream in=new GZIPInputStream(new BufferedInputStream(Files.newInputStream(path)),65536)){
				for(int size=in.read(buffer); size!=-1; size=in.read(buffer)){decoded.update(buffer,0,size); n+=size;}
			}
			require(n==bytes && DigestSuffix.fromDigest(decoded.digest()).equals(decodedPin),"BGZF9 decoded bytes differ from emitted rows");
		}
		final Path path; final BgzfOutputStream out; final MessageDigest digest;
		long bytes,members,residues; String decodedPin;
	}

	private static ArrayList<String[]> table(Path path,String hash,String header,int width) throws Exception{
		pin(path,hash); final ArrayList<String[]> rows=new ArrayList<String[]>();
		final ByteFile bf=ByteFile.makeByteFile1(path.toString(),false); final LineParser1 lp=new LineParser1('\t');
		try{
			if(header!=null){require(Arrays.equals(bf.nextLine(),header.getBytes(StandardCharsets.US_ASCII)),"Wrong metadata header: "+path);}
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lp.set(line); require(lp.terms()==width,"Wrong metadata width: "+path); final String[] r=new String[width];
				for(int i=0; i<width; i++){r[i]=lp.parseString(i); require(!r[i].isEmpty(),"Empty metadata field");} rows.add(r);
			}
		}finally{require(!bf.close(),"Metadata read failed: "+path);}
		pin(path,hash); return rows;
	}
	private static HashMap<String,String> kv(Path path,String hash) throws Exception{
		final HashMap<String,String> values=new HashMap<String,String>();
		for(String[] r : table(path,hash,null,2)){require(values.put(r[0],r[1])==null,"Duplicate metadata key");} return values;
	}
	private static long[] summary(Path path,String hash) throws Exception{
		final ArrayList<String[]> rows=table(path,hash,SUMMARY_HEADER,8); require(rows.size()==1,"Summary must have one row");
		final long[] n=new long[8]; for(int i=0; i<n.length; i++){n[i]=number(rows.get(0)[i]);}
		require(n[2]<=n[0] && n[3]==n[0]-n[2] && n[4]==Math.multiplyExact(n[0],200) && n[5]<=n[4] && n[6]<=n[0] && n[7]==n[0],"Invalid winner summary conservation"); return n;
	}
	private static HashMap<String,String> readPass(Path dir,String marker) throws Exception{
		final java.util.List<String> lines=Files.readAllLines(regular(dir.resolve("PASS")),StandardCharsets.US_ASCII);
		require(!lines.isEmpty() && lines.get(0).equals(marker),"Missing partition PASS"); final HashMap<String,String> pins=new HashMap<String,String>();
		for(int i=1; i<lines.size(); i++){
			final String[] r=lines.get(i).split("\t",-1); require(r.length==2 && r[0].matches("[A-Za-z_]+[.]tsv") && pins.put(r[0],r[1])==null,"Bad artifact pin row"); pin(dir.resolve(r[0]),r[1]);
		}
		return pins;
	}
	private static void publish(Path out,String marker,String[] files) throws Exception{
		final StringBuilder text=new StringBuilder(marker).append('\n');
		for(String file : files){text.append(file).append('\t').append(DigestSuffix.file(out.resolve(file).toString())).append('\n');}
		write(out.resolve("PASS"),text.toString());
	}
	private static void writeCounts(Path path,String[] ids,long[] counts) throws Exception{
		assert(ids.length==counts.length) : "Dense family identities and counts retain zero-member rows";
		final ByteBuilder text=new ByteBuilder(); text.append(FAMILY_HEADER).nl();
		for(int i=0; i<ids.length; i++){text.append(i).tab().append(ids[i]).tab().append(counts[i]).nl();}
		Files.write(path,Arrays.copyOf(text.array,text.length()));
	}
	private static String numbers(long[] values){final StringBuilder s=new StringBuilder(); for(long n : values){if(s.length()>0){s.append('\t');} s.append(n);} return s.toString();}
	private static void write(Path path,String text) throws Exception{Files.write(path,text.getBytes(StandardCharsets.US_ASCII));}
	private static Path regular(Path path) throws Exception{
		require(path.isAbsolute() && path.normalize().equals(path) && Files.isRegularFile(path,LinkOption.NOFOLLOW_LINKS) && path.toRealPath().equals(path),"Expected canonical regular file: "+path); return path;
	}
	private static void pin(Path path,String expected) throws Exception{DigestSuffix.requireSuffix(expected,"input pin"); require(DigestSuffix.file(regular(path).toString()).equals(expected),"Stored-byte pin mismatch: "+path);}
	private static long number(String text){final long n=Long.parseLong(text); require(n>=0 && Long.toString(n).equals(text),"Expected canonical nonnegative integer"); return n;}
	private static int integer(LineParser1 parser,int term){
		parser.setBounds(term); final byte[] line=parser.line(); final int start=parser.a(), end=parser.b();
		require(end>start && (end-start==1 || line[start]!='0'),"Expected canonical nonnegative assignment integer");
		int value=0;
		for(int i=start; i<end; i++){
			final int digit=line[i]-'0'; require(digit>=0 && digit<=9 && value<=(Integer.MAX_VALUE-digit)/10,"Malformed or overflowing assignment integer");
			value=value*10+digit;
		}
		return value;
	}
	private static String required(HashMap<String,String> o,String key){return HmmComparisonData.required(o,key);}
	private static boolean safe(String s){return s.matches("[A-Za-z0-9][A-Za-z0-9_.-]{0,127}");}
	private static boolean token(String s){if(s.isEmpty()){return false;} for(int i=0; i<s.length(); i++){if(s.charAt(i)<=32 || s.charAt(i)>=127){return false;}} return true;}
	private static void require(boolean ok,String message){if(!ok){throw new IllegalArgumentException(message);}}
	public static final String MEMBER_HEADER="family_rank\tfamily_id\tquery_id\tgenome_id\tsource_ordinal\tsequence";
	public static final String CHUNK_HEADER="batch_id\tbucket\tfile\tsha80\tmembers\tresidues";
	public static final String FAMILY_HEADER="roster_rank\tfamily_id\tassigned";
	private static final String BATCH_HEADER="batch_id\tshards\tshards_sha80\tresult\trecipe_sha80\tsummary_sha80";
	private static final String SHARD_HEADER="genome_id\tqueries\tresidues\tinput\tinput_sha80";
	private static final String SUMMARY_HEADER="queries\tresidues\tassigned\trejected\thbm_scored\thbm_nonpositive\ttied_winners\tblosum_scored";
	private static final String SUPPORT_HEADER="rank\tfamily_id\tconstruction_members\tassigned_members\tassigned_minus_construction";
	private HbmWinnerPartition(){}
}
