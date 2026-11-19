package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;

import ddl.DDLFamilyHoldout;
import ddl.DDLFamilySearch.StructuralInputs;
import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import prot.cluster.ClusterSketchCodec;
import prot.cluster.ClusteringInput;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.ListNum;

/** Converts final DDLClust member assignments into FamilyCalibration reference buckets.
 * Source FASTAs are streamed in declared DDLClust input order, avoiding a corpus-sized ID join.
 * @author Yelan
 */
public final class DdlClustCalibrationAdapter{

	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args, DdlClustCalibrationAdapter.class, false);
			try{runForTest(pp.args);}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	/** Throwing entry for fixtures and shell execution. */
	public static void runForTest(String[] args) throws Exception{
		Shared.setThreads(1); Shared.AMINO_IN=true; Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=Read.LOWER_CASE_TO_N=Read.DOT_DASH_X_TO_N=Read.FIX_HEADER=false;
		Read.JUNK_MODE=Read.IGNORE_JUNK;
		ReadWrite.ZIPLEVEL=9; ReadWrite.setZipThreads(1);
		ReadWrite.USE_BGZF=ReadWrite.PREFER_BGZIP=true; ReadWrite.USE_BGZIP=false;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final HashMap<String,String> o=options(args);
		final String verification=o.getOrDefault("verification", "structural");
		require(verification.equals("sha80") || verification.equals("structural"), "verification must be sha80 or structural");
		final StructuralInputs structural=verification.equals("structural") ? new StructuralInputs() : null;
		final int compressionThreads=Integer.parseInt(o.getOrDefault("compressionthreads","1"));
		require(compressionThreads>=1 && compressionThreads<=96,"compressionthreads must be1 through96");
		final boolean keepWork=Parse.parseBoolean(o.getOrDefault("keepwork", "f"));
		final boolean hasResidentRows=o.containsKey("maxresidentrows");
		final boolean hasResidentBytes=o.containsKey("maxresidentbytes");
		require(hasResidentRows==hasResidentBytes,
			"Resident storage requires maxresidentrows and maxresidentbytes together");
		require(!keepWork, "Resident storage has no work files to retain");
		final Long explicitMaxRows=hasResidentRows ? number(required(o,"maxresidentrows")) : null;
		final Long explicitMaxBytes=hasResidentBytes ? number(required(o,"maxresidentbytes")) : null;
		final MemberLengthRange memberRange=MemberLengthRange.from(o);
		final int tailCrop=Integer.parseInt(o.getOrDefault("tailcrop","0"));
		require(tailCrop==0 || tailCrop==100,"Adapter tailcrop must be0 or100");
		require(tailCrop==0 || memberRange.contains(tailCrop),"Cropped tail must fit the declared member range");
		final Path members=regular(Paths.get(required(o, "members")).toAbsolutePath().normalize());
		final String[] inputNames=ClusteringInput.inputFiles(required(o, "input"));
		final Path[] fasta=new Path[inputNames.length];
		for(int i=0; i<fasta.length; i++){
			fasta[i]=regular(Paths.get(inputNames[i]).toAbsolutePath().normalize());
		}
		final Path codec=regular(Paths.get(required(o, "codec")).toAbsolutePath().normalize());
		final Path out=Paths.get(required(o, "out")).toAbsolutePath().normalize();
		require(!Files.exists(out, LinkOption.NOFOLLOW_LINKS) && Files.isDirectory(out.getParent()), "Output must be fresh under an existing parent: "+out);
		final int buckets=Math.toIntExact(number(o.getOrDefault("buckets", "256")));
		require(buckets>0 && buckets<=256, "Bucket count out of range: "+buckets);
		final String membersPin=identity(members, structural);
		final String[] fastaPins=new String[fasta.length];
		for(int i=0; i<fasta.length; i++){fastaPins[i]=identity(fasta[i], structural);}
		final String codecPin=identity(codec, structural);
		final Family[] families=loadFamilies(codec, optional(o, "expectedfamilies"));
		final IntHashMap familyById=new IntHashMap(families.length*2+1);
		long expected=0;
		for(int i=0; i<families.length; i++){
			final Family family=families[i];
			require(familyById.get(family.stableId)<0, "Duplicate stable family ID: "+family.stableId);
			familyById.put(family.stableId, i);
			expected=Math.addExact(expected, family.assigned);
		}
		final Long expectedMembers=optional(o, "expectedmembers");
		require(expectedMembers==null || expectedMembers.longValue()==expected, "Expected member count differs from family roster");
		require(explicitMaxRows==null || explicitMaxRows.longValue()>0 && explicitMaxBytes.longValue()>0,
			"Resident storage limits must be positive");
		final long maxRows=explicitMaxRows==null ? expected : explicitMaxRows.longValue();
		final long maxBytes=explicitMaxBytes==null ? defaultResidentBytes() : explicitMaxBytes.longValue();
		final String budgetPolicy=explicitMaxBytes==null ? "auto_heap_available_advanced_minus_max_64MiB_or_25pct" : "explicit";
		final ResidentBuckets resident=new ResidentBuckets(families,buckets,expected,maxRows,maxBytes,compressionThreads,budgetPolicy);
		requireCompressors(compressionThreads);
		Files.createDirectory(out);
		final long[] routed=new long[buckets], routedResidues=new long[buckets];
		final long[] skipped=routeAssignments(members, fasta, families, familyById, resident, routed, routedResidues, buckets, memberRange, tailCrop);
		for(Family family : families){
			require(family.observed==family.assigned, "Family count mismatch for "+family.stableId+"/"+family.id+
				": observed="+family.observed+" expected="+family.assigned);
		}
		final Path[] bucketFiles=new Path[buckets];
		final long[] bucketMembers=new long[buckets], bucketResidues=new long[buckets];
		for(int bucket=0; bucket<buckets; bucket++){
			bucketFiles[bucket]=out.resolve(routeName("bucket", bucket));
			final long[] counts=resident.writeBucket(bucketFiles[bucket],bucket,routed[bucket],routedResidues[bucket]);
			bucketMembers[bucket]=counts[0]; bucketResidues[bucket]=counts[1];
		}
		require(identity(members, structural).equals(membersPin), "Member table changed while adapter was running");
		for(int i=0; i<fasta.length; i++){
			require(identity(fasta[i], structural).equals(fastaPins[i]), "Input FASTA changed while adapter was running: "+fasta[i]);
		}
		require(identity(codec, structural).equals(codecPin), "Cluster codec changed while adapter was running");
		writeOutputs(out, members, membersPin, fasta, fastaPins, codec, codecPin, families, bucketFiles,
			bucketMembers, bucketResidues, buckets, expected, skipped, memberRange, verification, structural, resident,
			tailCrop);
		System.err.println("DDLCLUST_CALIBRATION_ADAPTER_PASS families="+families.length+" members="+expected+" buckets="+buckets);
	}

	private static long[] routeAssignments(Path membersPath, Path[] fasta, Family[] families, IntHashMap familyById,
			ResidentBuckets resident, long[] routed, long[] residues, int buckets, MemberLengthRange memberRange, int tailCrop){
		final ByteFile table=ByteFile.makeByteFile1(membersPath.toString(), false);
		final ProteinCursor proteins=new ProteinCursor(fasta);
		final LineParser1 lp=new LineParser1('\t');
		final ByteBuilder row=new ByteBuilder();
		long assigned=0, unassigned=0, filtered=0;
		try{
			require(Arrays.equals(table.nextLine(), CLUSTERING_MEMBER_HEADER), "Wrong DDLClust member header");
			for(int ordinal=0; ; ordinal++){
				final Read read=proteins.next();
				if(read==null){
					require(table.nextLine()==null, "Member table contains more rows than the input FASTA");
					break;
				}
				final byte[] line=table.nextLine();
				require(line!=null, "Member table ended before the input FASTA");
				lp.set(line);
				require(lp.terms()==4 && lp.parseInt(0)==ordinal, "Member row is malformed or out of input order");
				final String id=firstToken(read.id);
				require(lp.termEquals(id, 1), "Member row input_id differs from FASTA ordinal "+ordinal+": "+id);
				final byte[] sequence=tailCrop==0 ? read.bases : croppedTail(read.bases,tailCrop,ordinal);
				final int stableId=lp.parseInt(2);
				if(lp.termEquals("assigned", 3)){
					require(stableId>=0, "Assigned row has a negative stable family ID");
					final int rank=familyById.get(stableId);
					require(rank>=0, "Assigned row references a stable family ID absent from the codec: "+stableId);
					final Family family=families[rank];
					validateAssigned(id, sequence, memberRange);
					family.observed++; family.residues=Math.addExact(family.residues, sequence.length);
					final int bucket=family.rank%buckets;
					row.clear().append(family.rank).tab().append(family.id).tab().append(id).tab();
					final int genomeStart=row.length();
					row.append(genome(id));
					final int genomeLength=row.length()-genomeStart;
					row.tab().append(ordinal).tab().append(sequence).nl();
					resident.add(bucket,row,genomeStart,genomeLength,ordinal,sequence.length);
					routed[bucket]++; residues[bucket]=Math.addExact(residues[bucket], sequence.length); assigned++;
				}else if(lp.termEquals("unassigned", 3)){
					require(stableId==-1, "Unassigned row has a family ID: "+stableId);
					unassigned++;
				}else if(lp.termEquals("filtered_length", 3)){
					require(stableId==-1, "Length-filtered row has a family ID: "+stableId);
					filtered++;
				}else{throw new IllegalArgumentException("Unknown DDLClust member status at ordinal "+ordinal);}
			}
		}finally{
			final boolean tableFailed=table.close();
			proteins.close();
			require(!tableFailed, "Member table read failed");
		}
		return new long[]{assigned, unassigned, filtered};
	}

	private static Family[] loadFamilies(Path path, Long expectedFamilies){
		final ClusterSketchCodec.Loaded loaded=ClusterSketchCodec.read(path.toString());
		require(loaded.header.families==loaded.families.length && loaded.families.length>0, "Codec family count mismatch or empty roster");
		require(expectedFamilies==null || expectedFamilies.longValue()==loaded.families.length, "Expected family count differs from codec");
		final Family[] families=new Family[loaded.families.length];
		for(int rank=0; rank<families.length; rank++){
			final ClusterSketchCodec.Entry entry=loaded.families[rank];
			require(safeFamily(entry.name), "Codec family name is not safe for FamilyCalibration: "+entry.name);
			families[rank]=new Family(rank, entry.familyId, entry.name, entry.members, entry.consensus);
		}
		return families;
	}

	private static void writeOutputs(Path out, Path members, String membersPin, Path[] fasta, String[] fastaPins, Path codec, String codecPin,
			Family[] families, Path[] buckets, long[] bucketMembers, long[] bucketResidues, int bucketCount, long assigned, long[] skipped,
			MemberLengthRange memberRange, String verification, StructuralInputs structural, ResidentBuckets resident,
			int tailCrop) throws Exception{
		final ByteBuilder line=new ByteBuilder();
		ByteStreamWriter writer=HmmComparisonData.writer(out.resolve("families.tsv").toString());
		try{
			writer.println(DDLFamilyHoldout.FAMILY_HEADER);
			for(Family family : families){writer.println(line.clear().append(family.rank).tab().append(family.id).tab().append(family.assigned));}
		}finally{HmmComparisonData.close(writer);}
		writer=HmmComparisonData.writer(out.resolve("family_map.tsv").toString());
		try{
			writer.println("stable_family_id\troster_rank\tfamily_id\tassigned\tresidues");
			for(Family family : families){
				writer.println(line.clear().append(family.stableId).tab().append(family.rank).tab()
					.append(family.id).tab().append(family.assigned).tab().append(family.residues));
			}
		}finally{HmmComparisonData.close(writer);}
		writer=HmmComparisonData.writer(out.resolve("consensus.faa").toString());
		try{
			for(Family family : families){
				line.clear().append('>').append(family.id).nl();
				appendConsensus(line, family.consensus).nl();
				writer.print(line);
			}
		}finally{HmmComparisonData.close(writer);}
		final String consensusPin=identity(out.resolve("consensus.faa"), structural);
		writer=HmmComparisonData.writer(out.resolve("models.tsv").toString());
		try{
			writer.println(HbmTailLibrary.MODELS_HEADER);
			for(Family family : families){
				writer.println(line.clear().append(family.rank).tab().append(family.id).tab().append(out.toString())
					.tab().append(family.assigned).tab().append(family.residues).tab()
					.append(".\t").append(consensusPin).append("\t.\t."));
			}
		}finally{HmmComparisonData.close(writer);}
		writer=HmmComparisonData.writer(out.resolve("chunks.tsv").toString());
		try{
			writer.println(DDLFamilyHoldout.CHUNK_HEADER);
			for(int bucket=0; bucket<bucketCount; bucket++){
				writer.println(line.clear().append("b00").tab().append(bucket).tab().append(buckets[bucket].toString())
					.tab().append(identity(buckets[bucket], structural)).tab().append(bucketMembers[bucket]).tab().append(bucketResidues[bucket]));
			}
		}finally{HmmComparisonData.close(writer);}
		writer=HmmComparisonData.writer(out.resolve("buckets.tsv").toString());
		try{
			writer.println("bucket\tmembers\tresidues\tfile");
			for(int bucket=0; bucket<bucketCount; bucket++){
				writer.println(line.clear().append(bucket).tab().append(bucketMembers[bucket]).tab()
					.append(bucketResidues[bucket]).tab().append(buckets[bucket].toString()));
			}
		}finally{HmmComparisonData.close(writer);}
		assert(fasta.length>0 && fasta.length==fastaPins.length) : "Input provenance requires one identity per ordered FASTA";
		final StringBuilder inputs=new StringBuilder("input\tpath\tsha80\nmembers\t").append(members).append('\t').append(membersPin).append('\n');
		for(int i=0; i<fasta.length; i++){
			inputs.append("fasta");
			if(fasta.length>1){inputs.append(i);}
			inputs.append('\t').append(fasta[i]).append('\t').append(fastaPins[i]).append('\n');
		}
		inputs.append("codec\t").append(codec).append('\t').append(codecPin).append('\n');
		write(out.resolve("inputs.sha80.tsv"), inputs.toString());
		write(out.resolve("summary.tsv"), "key\tvalue\nfamilies\t"+families.length+"\nassigned_members\t"+assigned+
			"\nunassigned_members\t"+skipped[1]+"\nfiltered_length_members\t"+skipped[2]+"\nbuckets\t"+bucketCount+"\n");
		write(out.resolve("recipe.tsv"), "key\tvalue\nschema\tddlclust_calibration_adapter_v1\nsource_order\tddlclust_input_index\n"+
			"family_map\tClusterSketchCodec stable numeric ID to display name\nbucket_order\tgenome_source_ordinal_query\nbatch_id\tb00\n"+memberRange.recipe("\t")+
			(tailCrop==0 ? "" : "tailcrop\t"+tailCrop+"\nterminal_stop\tstrip_one\nretained_residues\tstrict\n")+
			"verification\t"+verification+"\n"+(structural==null ? "" : "identity\tunmeasured\n")+
			resident.recipe());
		if(structural!=null){structural.verify();}
		write(out.resolve("PASS"), "DDLCLUST_CALIBRATION_ADAPTER_PASS\n");
	}

	/** Structural mode records stability, never an invented content identity. External
	 * guards bind the accepted inputs; all decoded conservation checks remain active. */
	private static String identity(Path file, StructuralInputs structural) throws Exception{
		if(structural==null){return DigestSuffix.file(file.toString());}
		structural.track(file);return "unmeasured";
	}

	private static HashMap<String,String> options(String[] args){
		final HashMap<String,String> out=new HashMap<String,String>();
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && pair[1]!=null && !pair[1].isEmpty(), "Expected flag=value: "+arg);
			final String key=pair[0].equalsIgnoreCase("in") ? "input" : pair[0].toLowerCase(Locale.ROOT);
			require(Arrays.asList("members", "input", "codec", "out", "buckets", "expectedfamilies",
				"expectedmembers", "keepwork", "minmemberaa", "maxmemberaa", "verification",
				"maxresidentrows", "maxresidentbytes", "tailcrop", "compressionthreads").contains(key),
				"Unknown adapter option: "+key);
			require(out.put(key, pair[1])==null, "Duplicate option: "+key);
		}
		return out;
	}

	/** Match ClusteringInput's biological suffix before strict residue validation.
	 * Validate every original row, including unassigned rows, without filtering ordinals.
	 * Disabled cropping bypasses this method and preserves historical bytes and failures.
	 */
	private static byte[] croppedTail(byte[] raw,int length,int ordinal){
		require(raw!=null,"Null residues at input ordinal "+ordinal);
		final int end=raw.length-(raw.length>0 && raw[raw.length-1]=='*' ? 1 : 0);
		require(end>=length,"Protein at input ordinal "+ordinal+" is below tailcrop="+length+"; cropping must not filter the population");
		final int start=end-length;
		for(int i=start;i<end;i++){
			require(VALID_ASSIGNED[raw[i]&255]!=0,"Unscorable retained residue at input ordinal "+ordinal);
		}
		return start==0 && end==raw.length ? raw : Arrays.copyOfRange(raw,start,end);
	}

	private static long defaultResidentBytes(){
		final long available=Shared.memAvailableAdvanced();
		final long headroom=Math.max(64L<<20, available>>2);
		require(available>headroom,
			"Default resident storage has no heap after reserving "+headroom+" bytes of headroom");
		return available-headroom;
	}

	private static void validateAssigned(String id, byte[] sequence, MemberLengthRange memberRange){
		require(id!=null && sequence!=null, "Assigned row lacks an ID or sequence");
		final String genome=genome(id);
		validateAssignedId(id, genome);
		if(!memberRange.contains(sequence.length)){
			throw new IllegalArgumentException("Assigned sequence length outside "+memberRange+
				": "+id+" length="+sequence.length);
		}
		for(byte b : sequence){
			require(VALID_ASSIGNED[b&255]!=0, "Assigned sequence contains an unscorable residue: "+id);
		}
	}

	private static void validateAssignedId(String id, String genome){
		require(safeToken(id) && safeToken(genome), "Unsafe query or genome ID");
	}

	private static Path regular(Path path){
		require(Files.isRegularFile(path), "Missing regular input file: "+path);
		return path;
	}

	private static void write(Path path, String text) throws Exception{
		require(!Files.exists(path, LinkOption.NOFOLLOW_LINKS), "Refusing existing output: "+path);
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII));
	}

	private static String genome(String id){
		final int pipe=id.indexOf('|');
		require(pipe>0, "Query ID lacks a genome prefix: "+id);
		return id.substring(0, pipe);
	}

	private static String firstToken(String id){
		require(id!=null && !id.isEmpty(), "Missing FASTA header");
		int end=0;
		while(end<id.length() && id.charAt(end)>32){end++;}
		return id.substring(0, end);
	}

	private static String routeName(String prefix, int route){
		return prefix+"_"+String.format(Locale.ROOT, "%03d", route)+(prefix.equals("bucket") ? ".tsv.gz" : ".tsv");
	}

	private static String required(HashMap<String,String> o, String key){
		final String value=o.get(key); require(value!=null, "Missing "+key+"="); return value;
	}
	private static Long optional(HashMap<String,String> o, String key){return o.containsKey(key) ? number(o.get(key)) : null;}
	private static long number(String value){
		final long n=Long.parseLong(value);
		require(n>=0 && Long.toString(n).equals(value), "Expected nonnegative canonical integer: "+value);
		return n;
	}
	private static boolean safeFamily(String s){return s!=null && s.matches("[A-Za-z0-9_.-]+");}
	private static boolean safeToken(String s){
		if(s==null || s.isEmpty()){return false;}
		for(int i=0; i<s.length(); i++){if(s.charAt(i)<=32 || s.charAt(i)>=127){return false;}}
		return true;
	}
	private static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}

	private static ByteBuilder appendConsensus(ByteBuilder out, byte[] consensus){
		for(byte residue : consensus){out.append(residue==Blosum62.X_CODE ? (byte)'X' : AminoAcid.numberToAcid[residue]);}
		return out;
	}

	private static byte[] validAssigned(){
		final byte[] valid=new byte[256], single=new byte[1];
		for(int c=0; c<valid.length; c++){
			single[0]=(byte)c;
			try{Blosum62.encode(single, "assigned_validity"); valid[c]=1;}
			catch(RuntimeException ignored){}
		}
		return valid;
	}

	private static final byte[] CLUSTERING_MEMBER_HEADER="input_index\tinput_id\tfamily_id\tstatus".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] VALID_ASSIGNED=validAssigned();

	private static final class Family{
		Family(int rank_, int stableId_, String id_, int assigned_, byte[] consensus_){
			rank=rank_; stableId=stableId_; id=id_; assigned=assigned_; consensus=consensus_;
		}
		final int rank, stableId, assigned;
		final String id;
		final byte[] consensus;
		int observed;
		long residues;
	}

	/** Compact retained rows; genome bytes share the output row rather than a second String. */
	private static final class ResidentRow implements Comparable<ResidentRow>{
		ResidentRow(byte[] line_,int start_,int length_,long ordinal_,int residues_){
			line=line_;start=start_;length=length_;ordinal=ordinal_;residues=residues_;
		}
		/** IDs were validated as ASCII, so byte order matches the legacy String order. */
		int compareGenome(ResidentRow other){
			final int common=Math.min(length,other.length);
			for(int i=0;i<common;i++){
				final int difference=line[start+i]-other.line[other.start+i];
				if(difference!=0){return difference;}
			}
			return Integer.compare(length,other.length);
		}
		/** Source ordinals are unique; no query-ID tie breaker can affect a valid row. */
		@Override public int compareTo(ResidentRow other){
			final int order=compareGenome(other);
			return order!=0 ? order : Long.compare(ordinal,other.ordinal);
		}
		final byte[] line;
		final int start,length,residues;
		final long ordinal;
	}

	/** Bounded storage with no spill path. Accounting reserves 8-byte references,
	 * 32-byte array headers, 64-byte row objects, aligned payloads and a full largest-
	 * bucket sort workspace. This bounds retained row storage, not total JVM heap:
	 * codec/roster, stream buffers and the VM need separate headroom. */
	private static final class ResidentBuckets{
		/** Reserve exact roster-sized pointer arrays and sort space before accepting rows. */
		ResidentBuckets(Family[] families,int buckets,long expected,long maxRows_,long maxBytes_,int compressionThreads_,String budgetPolicy_){
			maxRows=maxRows_;maxBytes=maxBytes_;compressionThreads=compressionThreads_;budgetPolicy=budgetPolicy_;
			require(maxRows>=0 && expected<=maxRows,"Assigned roster exceeds maxresidentrows: "+expected+" > "+maxRows);
			require(maxBytes>0,"maxresidentbytes must be positive");
			final long[] sizes=new long[buckets];
			for(Family family:families){sizes[family.rank%buckets]=Math.addExact(sizes[family.rank%buckets],family.assigned);}
			long largest=0;
			for(long size:sizes){require(size<=Integer.MAX_VALUE-8,"Resident bucket exceeds array capacity");largest=Math.max(largest,size);}
			charge(arrayBytes(buckets, expected, largest));
			rows=new ResidentRow[buckets][];used=new int[buckets];rawBytes=new long[buckets];
			for(int i=0;i<buckets;i++){
				rows[i]=new ResidentRow[(int)sizes[i]];
				rawBytes[i]=DDLFamilyHoldout.MEMBER_HEADER.length()+1L;
			}
		}
		/** Charge before retaining any payload so over-budget data never silently spill. */
		void add(int bucket,ByteBuilder line,int start,int length,long ordinal,int residues){
			require(used[bucket]<rows[bucket].length,"Assigned bucket exceeds codec member count: "+bucket);
			require(count<maxRows,"Assigned rows exceed maxresidentrows");
			charge(96L+((line.length()+7L)&~7L));
			rows[bucket][used[bucket]++]=new ResidentRow(line.toBytes(),start,length,ordinal,residues);
			count++;rawBytes[bucket]=Math.addExact(rawBytes[bucket],line.length());
		}
		/** Account conservatively before allocation and reject the first excess. */
		void charge(long bytes){
			accountedBytes=Math.addExact(accountedBytes,bytes);
			require(accountedBytes<=maxBytes,"Resident storage exceeds maxresidentbytes: "+accountedBytes+" > "+maxBytes);
		}
		static long arrayBytes(Family[] families,int buckets){
			final long[] sizes=new long[buckets];
			long expected=0, largest=0;
			for(Family family:families){
				sizes[family.rank%buckets]=Math.addExact(sizes[family.rank%buckets], family.assigned);
				expected=Math.addExact(expected, family.assigned);
			}
			for(long size:sizes){largest=Math.max(largest, size);}
			return arrayBytes(buckets, expected, largest);
		}
		static long arrayBytes(int buckets,long expected,long largest){
			return Math.addExact(1024L+128L*buckets,Math.multiplyExact(8L,Math.addExact(expected, largest)));
		}
		/** Sort one exact-sized bucket and write its final compressed representation. */
		long[] writeBucket(Path path,int bucket,long expectedRows,long expectedResidues){
			final ResidentRow[] data=rows[bucket];
			require(used[bucket]==data.length && data.length==expectedRows,"Resident bucket member conservation failed: "+bucket);
			Arrays.sort(data);
			final ByteStreamWriter writer=compressedWriter(path,rawBytes[bucket],compressionThreads);
			ResidentRow previous=null;
			long residues=0,bytes=DDLFamilyHoldout.MEMBER_HEADER.length()+1L;
			try{
				writer.println(DDLFamilyHoldout.MEMBER_HEADER);
				for(ResidentRow row:data){
					if(previous!=null){
						final int order=previous.compareGenome(row);
						require(order<0 || order==0 && previous.ordinal<row.ordinal,"Duplicate or reordered source ordinal in resident bucket");
					}
					writer.print(row.line);previous=row;
					residues=Math.addExact(residues,row.residues);bytes=Math.addExact(bytes,row.line.length);
				}
			}finally{HmmComparisonData.close(writer);}
			require(residues==expectedResidues && bytes==rawBytes[bucket],"Resident bucket residue/byte conservation failed: "+bucket);
			rows[bucket]=null;
			return new long[]{data.length,residues};
		}
		/** Record storage limits and compression policy without claiming a measured heap minimum. */
		String recipe(){
			return "storage\tbounded_resident\nmaxresidentrows\t"+maxRows+"\nmaxresidentbytes\t"+maxBytes+
				"\nresident_budget\t"+budgetPolicy+
				"\naccounted_resident_bytes\t"+accountedBytes+"\ncompression_raw_byte_boundary\t200000\n"+
				"compression_above_boundary\tbgzip9\ncompression_at_or_below_boundary\tpigz11\ncompression_threads\t"+compressionThreads+"\n";
		}
		final int compressionThreads;
		final ResidentRow[][] rows;
		final int[] used;
		final long[] rawBytes;
		final long maxRows,maxBytes;
		final String budgetPolicy;
		long count,accountedBytes;
	}

	/** Require exact codecs before output creation; native BGZF cannot provide level 9. */
	private static void requireCompressors(int compressionThreads){
		require(Data.BGZIP() && Data.BGZIP_VERSION_levelFlag,"Structural output requires external bgzip with level support");
		require(compressionThreads==1 || Data.BGZIP_VERSION_threadsFlag,"Multiple compressionthreads require threaded bgzip");
		require(Data.PIGZ() && Data.PIGZ_VERSION_23plus,"Structural output requires external pigz with level 11 support");
	}

	/** Open one bucket stream, restoring codec and thread globals before starting its writer.
	 * Public for the compression/restoration fixture; callers must serialize stream construction. */
	public static ByteStreamWriter compressedWriter(Path path,long rawBytes,int compressionThreads){
		require(compressionThreads>=1 && compressionThreads<=96,"compressionthreads must be1 through96");
		require(rawBytes<=200000 || compressionThreads==1 || Data.BGZIP() && Data.BGZIP_VERSION_threadsFlag,
			"Multiple compressionthreads require threaded bgzip");
		require(!Files.exists(path,LinkOption.NOFOLLOW_LINKS),"Refusing existing compressed bucket: "+path);
		final int level=ReadWrite.ZIPLEVEL;
		final int threads=Shared.threads(),buffers=Shared.numBuffers(),zipThreads=ReadWrite.MAX_ZIP_THREADS();
		final float multiplier=ReadWrite.ZIP_THREAD_MULT();
		final boolean setThreads=ReadWrite.SET_ZIP_THREADS,setMultiplier=ReadWrite.SET_ZIP_THREAD_MULT;
		final boolean bgzf=ReadWrite.USE_BGZF,bgzip=ReadWrite.USE_BGZIP,pigz=ReadWrite.USE_PIGZ,
			forceBgzip=ReadWrite.FORCE_BGZIP,forcePigz=ReadWrite.FORCE_PIGZ,
			prefer=ReadWrite.PREFER_BGZIP,nativeBgzf=ReadWrite.ALLOW_NATIVE_BGZF;
		final ByteStreamWriter writer;
		try{
			Shared.setThreads(compressionThreads);ReadWrite.setZipThreads(compressionThreads);ReadWrite.setZipThreadMult(1f);
			final boolean large=rawBytes>200000;
			ReadWrite.ZIPLEVEL=large ? 9 : 11;
			ReadWrite.USE_BGZF=ReadWrite.USE_BGZIP=ReadWrite.FORCE_BGZIP=ReadWrite.PREFER_BGZIP=large;
			ReadWrite.USE_PIGZ=ReadWrite.FORCE_PIGZ=!large;ReadWrite.ALLOW_NATIVE_BGZF=false;
			writer=new ByteStreamWriter(path.toString(),false,false,true);
		}finally{
			Shared.setThreads(threads);Shared.setBuffers(buffers);
			ReadWrite.setZipThreads(zipThreads);ReadWrite.SET_ZIP_THREADS=setThreads;
			ReadWrite.setZipThreadMult(multiplier);ReadWrite.SET_ZIP_THREAD_MULT=setMultiplier;
			ReadWrite.ZIPLEVEL=level;ReadWrite.USE_BGZF=bgzf;ReadWrite.USE_BGZIP=bgzip;ReadWrite.USE_PIGZ=pigz;
			ReadWrite.FORCE_BGZIP=forceBgzip;ReadWrite.FORCE_PIGZ=forcePigz;ReadWrite.PREFER_BGZIP=prefer;
			ReadWrite.ALLOW_NATIVE_BGZF=nativeBgzf;
		}
		writer.start();return writer;
	}

	/** Holds one source stream at a time; file boundaries never reset the caller's ordinal. */
	private static final class ProteinCursor{
		ProteinCursor(Path[] paths_){
			require(paths_!=null && paths_.length>0, "Protein cursor requires an ordered nonempty input list");
			paths=paths_;openNext();
		}
		/** Starts the next declared source only after the prior stream has closed successfully. */
		private void openNext(){
			assert(stream==null && fileIndex<paths.length) : "Close the preceding FASTA before opening the next ordered source";
			final FileFormat format=FileFormat.testInput(paths[fileIndex++].toString(), FileFormat.FASTA, null, true, true);
			stream=StreamerFactory.makeStreamer(format, null, true, -1, false, true, 1);
			stream.start();
			batch=null;index=0;
		}
		Read next(){
			while(batch==null || index>=batch.size()){
				if(stream==null){return null;}
				final ListNum<Read> next=stream.nextList();
				if(next==null){
					close();
					if(fileIndex==paths.length){return null;}
					openNext();continue;
				}
				batch=next.list; index=0;
			}
			final Read read=batch.get(index++);
			require(read!=null && read.id!=null && read.bases!=null && read.mate==null, "Expected an unpaired protein FASTA record");
			return read;
		}
		void close(){
			if(stream==null){return;}
			final Streamer current=stream;stream=null;
			if(ReadWrite.closeStream(current)){throw new IllegalStateException("Protein FASTA input failed: "+paths[fileIndex-1]);}
		}
		final Path[] paths;
		Streamer stream;
		ArrayList<Read> batch;
		int index,fileIndex;
	}
}
