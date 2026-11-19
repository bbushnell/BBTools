package prot;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Objects;

import ddl.DDLFamilySearch;
import ddl.DDLKeyTieBreak;
import ddl.DDLKeyTieBreak.Policy;
import ddl.DDLKeyTieBreak.PolicyContext;
import ddl.DDLPentamerBlacklist;
import ddl.FamilySketchBuilder;
import ddl.ProteinBlacklistPolicy;
import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import shared.Shared;
import structures.IntHashMap;
import structures.ListNum;

/** Full-view raw counts or binding-only validation for resident clustering stages.
 * Reads only explicitly declared inputs, with no count sidecar or per-worker copy.
 * @author Yaoyao
 */
public final class FamilyKeyPolicyInputs{

	/** Build one immutable policy context; no sequence reads for non-raw policies. */
	public static PolicyContext load(Policy policy,Path rawManifest,long expectedRawTails,
			int rawMinLength,int rawMaxLength,int rawTailCrop,int width,ProteinKmerGeometry geometry,
			long sketchSeed,ProteinPositionMode positionMode,int partitions,int overlap,
			DDLPentamerBlacklist blacklist,DDLFamilySearch.StructuralInputs tracked,PrintStream log) throws Exception{
		return load(policy,rawManifest,expectedRawTails,rawMinLength,rawMaxLength,rawTailCrop,width,geometry,
			sketchSeed,positionMode,partitions,overlap,blacklist,tracked,log,null);
	}

	/** A resident context is checked against declared metadata, never recounted or reranked. */
	public static PolicyContext load(Policy policy,Path rawManifest,long expectedRawTails,
			int rawMinLength,int rawMaxLength,int rawTailCrop,int width,ProteinKmerGeometry geometry,
			long sketchSeed,ProteinPositionMode positionMode,int partitions,int overlap,
			DDLPentamerBlacklist blacklist,DDLFamilySearch.StructuralInputs tracked,PrintStream log,
			PolicyContext resident) throws Exception{
		return load(policy,rawManifest,expectedRawTails,rawMinLength,rawMaxLength,rawTailCrop,width,geometry,
			sketchSeed,positionMode,partitions,overlap,blacklist,tracked,log,resident,0);
	}

	/** Bonus metadata is family-only and must agree across every resident stage. */
	public static PolicyContext load(Policy policy,Path rawManifest,long expectedRawTails,
			int rawMinLength,int rawMaxLength,int rawTailCrop,int width,ProteinKmerGeometry geometry,
			long sketchSeed,ProteinPositionMode positionMode,int partitions,int overlap,
			DDLPentamerBlacklist blacklist,DDLFamilySearch.StructuralInputs tracked,PrintStream log,
			PolicyContext resident,long consensusMemberBonus) throws Exception{
		require(consensusMemberBonus>=0 && consensusMemberBonus<=DDLKeyTieBreak.MAX_VOTES,
			"Consensus member bonus must fit a nonnegative40-bit selection score");
		require(policy!=null,"Missing key policy");
		final DDLPentamerBlacklist selected=ProteinBlacklistPolicy.resolve(geometry,blacklist);
		final String pin=selected==null ? null : selected.sha80();
		if(resident!=null){
			resident.requireGeometry(width,geometry,sketchSeed,positionMode,partitions,overlap,pin);
			require(resident.policy==policy && resident.consensusMemberBonus==consensusMemberBonus &&
				resident.encodingVersion.equals(DDLKeyTieBreak.encodingVersion(policy)),
				"Resident policy/encoding differs from the requested recipe");
		}
		if(policy!=Policy.RAW_TAIL_RARITY){
			return resident!=null ? resident : new PolicyContext(policy,null,width,geometry,sketchSeed,positionMode,
				partitions,overlap,pin,null,consensusMemberBonus);
		}
		require(expectedRawTails>0 && rawMinLength>0 && rawMaxLength>=rawMinLength && rawTailCrop>=0 &&
			(rawTailCrop==0 || rawTailCrop>=rawMinLength && rawTailCrop<=rawMaxLength),
			"Raw counting requires a complete population and explicit compatible length/crop settings");
		final Manifest manifest=readManifest(rawManifest,tracked);
		final String binding="raw_input_v1\nexpected="+expectedRawTails+"\nmin="+rawMinLength+"\nmax="+rawMaxLength+
			"\ntailcrop="+rawTailCrop+"\n"+manifest.rows;
		if(resident!=null){
			require(Objects.equals(resident.rawInputBinding,binding),"Resident full-view raw input/crop binding differs");
			tracked.verify();return resident;
		}
		final long start=System.nanoTime();
		final FamilySketchBuilder kernel=new FamilySketchBuilder(width,geometry,sketchSeed,positionMode,partitions,overlap,selected);
		final DDLKeyTieBreak.RawTailCounter counter=new DDLKeyTieBreak.RawTailCounter(width,expectedRawTails);
		final RawScratch scratch=new RawScratch();
		final char[] sketch=new char[width];long processed=0;
		try(InputSettings settings=new InputSettings()){
			for(Path input:manifest.inputs){
				final FileFormat format=FileFormat.testInput(input.toString(),FileFormat.FASTA,null,true,true);
				require(format.fasta(),"Raw clustering view requires FASTA: "+input);
				final Streamer stream=StreamerFactory.makeStreamer(format,null,true,-1,false,true,1);stream.start();
				try{
					for(ListNum<Read> batch=stream.nextList();batch!=null;batch=stream.nextList()){
						for(Read read:batch.list){
							require(read!=null && read.id!=null && read.bases!=null && read.mate==null,"Expected an unpaired raw protein record");
							require(!read.id.isEmpty() && read.id.charAt(0)>32,"Raw record requires a nonempty first header token");
							final byte[] encoded=scratch.encode(read.bases,rawTailCrop,rawMinLength,rawMaxLength,processed);
							kernel.sketchEncoded(encoded,sketch);counter.add(sketch);processed++;
						}
					}
				}finally{require(!ReadWrite.closeStream(stream),"Raw protein input failed: "+input);}
				tracked.verify();
				if(log!=null){log.println("Raw full-view count: "+input+"; cumulative="+processed);}
			}
		}
		final long counted=System.nanoTime();
		final DDLKeyTieBreak.RawTailCounts counts=counter.freeze();
		tracked.verify();
		if(log!=null){
			log.println("Raw full-view count seconds="+((counted-start)*1e-9)+"; rank seconds="+((System.nanoTime()-counted)*1e-9)+
				"; tails="+counts.tails()+"; postings="+counts.postings()+"; max_count="+counts.maximumCount());
		}
		return new PolicyContext(policy,counts,width,geometry,sketchSeed,positionMode,partitions,overlap,pin,binding,consensusMemberBonus);
	}

	/** Binding-only resident validation; returns no replacement context or sequence data. */
	public static void validateResident(PolicyContext resident,Policy policy,Path rawManifest,long expectedRawTails,
			int rawMinLength,int rawMaxLength,int rawTailCrop,int width,ProteinKmerGeometry geometry,
			long sketchSeed,ProteinPositionMode positionMode,int partitions,int overlap,
			DDLPentamerBlacklist blacklist,DDLFamilySearch.StructuralInputs tracked) throws Exception{
		require(resident!=null,"Resident validation requires the original context");
		require(load(policy,rawManifest,expectedRawTails,rawMinLength,rawMaxLength,rawTailCrop,width,geometry,
			sketchSeed,positionMode,partitions,overlap,blacklist,tracked,null,resident)==resident,
			"Resident validation cannot replace immutable counts");
	}

	/** Explicit bonus counterpart; old resident callers retain W0. */
	public static void validateResident(PolicyContext resident,Policy policy,Path rawManifest,long expectedRawTails,
			int rawMinLength,int rawMaxLength,int rawTailCrop,int width,ProteinKmerGeometry geometry,
			long sketchSeed,ProteinPositionMode positionMode,int partitions,int overlap,
			DDLPentamerBlacklist blacklist,DDLFamilySearch.StructuralInputs tracked,long consensusMemberBonus) throws Exception{
		require(resident!=null,"Resident validation requires the original context");
		require(load(policy,rawManifest,expectedRawTails,rawMinLength,rawMaxLength,rawTailCrop,width,geometry,
			sketchSeed,positionMode,partitions,overlap,blacklist,tracked,null,resident,consensusMemberBonus)==resident,
			"Resident validation cannot replace immutable counts or consensus bonus");
	}

	/** Metadata-only accessor for output alias protection; ordinary/raw dual-input roles are valid. */
	public static Path[] declaredInputs(Path manifest,DDLFamilySearch.StructuralInputs tracked) throws Exception{
		final Manifest declared=readManifest(manifest,tracked);
		final Path[] paths=new Path[declared.inputs.size()+1];paths[0]=manifest.toAbsolutePath().normalize();
		for(int i=0;i<declared.inputs.size();i++){paths[i+1]=declared.inputs.get(i);}
		return paths;
	}

	/** Parse exactly the declared ordered file set and retain its structural bindings. */
	private static Manifest readManifest(Path manifest,DDLFamilySearch.StructuralInputs tracked) throws Exception{
		require(manifest!=null && tracked!=null,"Raw counting requires its manifest and lifetime structural tracker");
		manifest=tracked.track(manifest);
		final Manifest result=new Manifest();final HashSet<Path> unique=new HashSet<Path>();
		final StringBuilder rows=new StringBuilder("order\tpath\tbytes\tmtime_seconds\n");
		final ByteFile reader=ByteFile.makeByteFile(manifest.toString(),false);
		try{
			final byte[] header=reader.nextLine();
			require(header!=null && new String(header,StandardCharsets.US_ASCII).equals("order\tpath\tbytes\tmtime_seconds"),
				"Exact ordered raw-input manifest schema required");
			for(byte[] line=reader.nextLine();line!=null;line=reader.nextLine()){
				final String[] f=new String(line,StandardCharsets.UTF_8).split("\t",-1);
				require(f.length==4 && Integer.parseInt(f[0])==result.inputs.size(),"Require dense explicit raw-input order");
				final Path input=Paths.get(f[1]);
				require(input.isAbsolute() && input.normalize().equals(input) && unique.add(input),"Require distinct normalized absolute raw paths");
				tracked.track(input);final long bytes=Long.parseLong(f[2]),mtime=Long.parseLong(f[3]);
				require(bytes>=0 && Files.size(input)==bytes && Files.getLastModifiedTime(input).toMillis()/1000==mtime,
					"Raw input differs from declared metadata: "+input);
				rows.append(result.inputs.size()).append('\t').append(input).append('\t').append(bytes).append('\t').append(mtime).append('\n');
				result.inputs.add(input);
			}
		}finally{require(!reader.close(),"Raw input manifest read failed");}
		require(!result.inputs.isEmpty(),"Raw population must declare at least one complete input");
		tracked.verify();result.rows=rows.toString();return result;
	}

	/** Same suffix/stop/strict-encoding contract, with one reusable array per observed retained length. */
	private static final class RawScratch{
		final IntHashMap slots=new IntHashMap(64);
		final ArrayList<byte[]> buffers=new ArrayList<byte[]>();
		/** Crop first, then encode every retained residue; scratch is consumed before the next call. */
		byte[] encode(byte[] raw,int crop,int min,int max,long ordinal){
			final int end=raw.length-(raw.length>0 && raw[raw.length-1]=='*' ? 1 : 0);
			final int length=crop==0 ? end : crop,start=end-length;
			if(start<0 || length<min || length>max){
				throw new IllegalArgumentException("Raw view contains an out-of-contract record at ordinal "+ordinal+"; no silent filtering");
			}
			int slot=slots.get(length);
			if(slot<0){slot=buffers.size();slots.put(length,slot);buffers.add(new byte[length]);}
			final byte[] encoded=buffers.get(slot);
			for(int i=0;i<length;i++){
				final byte value=ENCODING[raw[start+i]&255];
				if(value<0){throw new IllegalArgumentException("Illegal retained raw residue at record "+ordinal+" offset "+i);}
				encoded[i]=value;
			}
			return encoded;
		}
	}
	/** Build strict lookup once from the canonical encoder, without per-record boxing or arrays. */
	private static byte[] encoding(){
		final byte[] table=new byte[256];Arrays.fill(table,(byte)-1);
		for(int c=0;c<128;c++){
			try{table[c]=Blosum62.encode(new byte[]{(byte)c},"raw_encoding_table")[0];}catch(RuntimeException unsupported){}
		}
		return table;
	}
	private static final byte[] ENCODING=encoding();
	/** Scoped setup before any worker/read phase; restore caller globals after all raw streams close. */
	private static final class InputSettings implements AutoCloseable{
		final boolean amino=Shared.AMINO_IN,trim=Shared.TRIM_READ_DESCRIPTION;
		final boolean upper=Read.TO_UPPER_CASE,lower=Read.LOWER_CASE_TO_N,dot=Read.DOT_DASH_X_TO_N,header=Read.FIX_HEADER;
		final int junk=Read.JUNK_MODE;
		InputSettings(){
			Shared.AMINO_IN=true;Shared.TRIM_READ_DESCRIPTION=false;
			Read.TO_UPPER_CASE=Read.LOWER_CASE_TO_N=Read.DOT_DASH_X_TO_N=Read.FIX_HEADER=false;
			Read.JUNK_MODE=Read.IGNORE_JUNK;
		}
		@Override public void close(){
			Shared.AMINO_IN=amino;Shared.TRIM_READ_DESCRIPTION=trim;
			Read.TO_UPPER_CASE=upper;Read.LOWER_CASE_TO_N=lower;Read.DOT_DASH_X_TO_N=dot;Read.FIX_HEADER=header;
			Read.JUNK_MODE=junk;
		}
	}
	private static final class Manifest{
		final ArrayList<Path> inputs=new ArrayList<Path>();
		String rows;
	}
	private static void require(boolean ok,String message){if(!ok){throw new IllegalArgumentException(message);}}
	private FamilyKeyPolicyInputs(){}
}
