package ddl;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;

import dna.AminoAcid;
import cardinality.DynamicDemiLog32;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.IntHashMap;
import map.LongHashMap;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Primitive, reusable kernels for the family-profile feasibility experiment.
 * The score encoding matches DynamicDemiLog.hashAndStore at exponent=6.
 * No identity estimator is applied to the modal profiles.
 * @author Sayu
 */
public final class DDLFamilyProfiles{

	private DDLFamilyProfiles(){}

	/** Build requested widths from every member of each selected manifest family.
	 * Shards select manifest ROW numbers, not family identifiers or active_index.
	 */
	public static void main(String[] args) throws Exception{
		args=new PreParser(args, DDLFamilyProfiles.class, false).args;
		String manifest=null, pin=null, namespace=null, out=null, sketch=SKETCH_PLAIN, aaBits=null, positionModeFlag=null;
		Boolean positionalFlag=null;
		boolean positions=false;
		int k=5, shard=0, shards=1, maxFamilies=Integer.MAX_VALUE, registerBits=16, exponent=6;
		int partitions=POSITIONAL_QUADRANTS, overlap=POSITIONAL_OVERLAP;
		int[] widths={128,256,512};
		long seed=0;
		final HashSet<String> seen=new HashSet<String>();
		for(String arg : args){
			final String[] parts=Parse.splitOnFirst(arg, '=');
			if(parts.length!=2){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String key=parts[0].toLowerCase(Locale.ROOT), value=parts[1];
			if(!seen.add(key)){throw new IllegalArgumentException("Duplicate flag: "+key);}
			if(key.equals("manifest")){manifest=value;}
			else if(key.equals("manifestsha80")){pin=value;}
			else if(key.equals("namespace")){namespace=value;}
			else if(key.equals("out")){out=value;}
			else if(key.equals("k")){k=Integer.parseInt(value);}
			else if(key.equals("aabits") || key.equals("aapartition")){aaBits=value;}
			else if(key.equals("widths")){widths=parseWidths(value);}
			else if(key.equals("seed")){seed=Long.parseLong(value);}
			else if(key.equals("registerbits")){registerBits=Integer.parseInt(value);}
			else if(key.equals("exponent")){exponent=Integer.parseInt(value);}
			else if(key.equals("sketch")){sketch=value.toLowerCase(Locale.ROOT);}
			else if(key.equals("positional")){positionalFlag=Parse.parseBoolean(value);}
			else if(key.equals("positionmode")){positionModeFlag=value;}
			else if(key.equals("positions")){positions=Parse.parseBoolean(value);}
			else if(key.equals("partitions")){partitions=Integer.parseInt(value);}
			else if(key.equals("overlap")){overlap=Integer.parseInt(value);}
			else if(key.equals("shard")){shard=Integer.parseInt(value);}
			else if(key.equals("shards")){shards=Integer.parseInt(value);}
			else if(key.equals("maxfamilies")){maxFamilies=Integer.parseInt(value);}
			else{throw new IllegalArgumentException("Unknown flag: "+key);}
		}
		if(seen.contains("aabits") && seen.contains("aapartition")){
			throw new IllegalArgumentException("Use either aabits or aapartition, not both");
		}
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(k, aaBits);
		final DDLPentamerBlacklist official=ProteinBlacklistPolicy.resolve(geometry,null);
		if(manifest==null || pin==null || !pin.matches("[0-9a-f]{20}") || namespace==null ||
			!namespace.matches("[A-Za-z0-9_-]+") || out==null || k<1 || k>12 ||
			shards<1 || shard<0 || shard>=shards || maxFamilies<1){
			throw new IllegalArgumentException("Require manifest, manifestsha80, namespace, fresh out; k=1..12, shard<shards");
		}
		final ProteinPositionMode positionMode=ProteinPositionMode.resolve(positionModeFlag,
			seen.contains("sketch") ? sketch : null,positionalFlag,ProteinPositionMode.PLAIN);
		sketch=positionMode.sketch;
		validateLayout(registerBits,exponent);
		validatePositionalLayout(partitions,overlap,geometry);
		if(positionMode.positional()){validatePositionalLayout(partitions,overlap,geometry,positionMode);}
		if(!positionMode.positional() && (seen.contains("partitions") || seen.contains("overlap"))){
			throw new IllegalArgumentException("partitions/overlap require positional sketches");
		}
		if(registerBits==32 && !geometry.isExact()){
			throw new IllegalArgumentException("Irregular amino-acid partitions are currently 16-bit only");
		}
		if(positionMode.positional() && (registerBits!=16 || seed!=0 ||
				widths.length!=1 || widths[0]!=POSITIONAL_BUCKETS)){
			throw new IllegalArgumentException("Positional sketch v1 requires registerbits=16 seed=0 widths=1024");
		}
		if(positions && (!positionMode.positional() || registerBits!=16)){
			throw new IllegalArgumentException("Position payloads require positional 16-bit profiles");
		}
		final Path manifestPath=Paths.get(manifest).toAbsolutePath();
		final byte[] hashBuffer=new byte[65536];
		requirePin(manifestPath, pin, hashBuffer);
		final Path output=Paths.get(out).toAbsolutePath();
		Files.createDirectory(output);
		Shared.setThreads(1);
		Shared.AMINO_IN=true;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final long started=System.nanoTime();
		final ByteFile table=ByteFile.makeByteFile(manifestPath.toString(), false);
		final byte[] header=table.nextLine();
		if(header==null){throw new IllegalArgumentException("Empty manifest");}
		final String[] columns=new String(header, StandardCharsets.US_ASCII).split("\t", -1);
		final boolean tails=Arrays.asList(columns).contains("directory");
		final int idColumn=column(columns,"family_id");
		final int pathColumn=column(columns,tails ? "directory" : "fasta_file");
		final int countColumn=column(columns,tails ? "members" : "assigned_count");
		final int pinColumn=column(columns,tails ? "retained_sha80" : "fasta_sha80");
		final int residuesColumn=tails ? column(columns,"residues") : -1;
		final ByteStreamWriter profiles=new ByteStreamWriter(output.resolve("profiles.tsv.gz").toString(),false,false,false);
		profiles.start();
		profiles.print(registerBits==16 ? (positions ? HEADER16_POSITIONS : HEADER16)+"\n" : HEADER32+"\n");
		final LineParser1 lp=new LineParser1('\t');
		final HashSet<String> familyIds=new HashSet<String>();
		final ByteBuilder lineOut=new ByteBuilder();
		int row=0, families=0;
		long members=0, residues=0, words=0;
		boolean complete=false;
		try{
			for(byte[] line=table.nextLine(); line!=null; line=table.nextLine()){
				lp.set(line);
				if(lp.terms()!=columns.length){throw new IllegalArgumentException("Manifest field count at row "+row);}
				final String id=lp.parseString(idColumn);
				if(id.length()==0 || !familyIds.add(id)){throw new IllegalArgumentException("Empty/duplicate family: "+id);}
				final int index=row++;
				if(index%shards!=shard || families>=maxFamilies){continue;}
				final int expected=lp.parseInt(countColumn);
				if(expected<1){throw new IllegalArgumentException("Nonpositive member count: "+id);}
				Path fasta=Paths.get(lp.parseString(pathColumn));
				if(tails){fasta=fasta.resolve("retained.faa");}
				else if(!fasta.isAbsolute()){fasta=manifestPath.getParent().resolve(fasta);}
				final String fastaPin=lp.parseString(pinColumn);
				if(!fastaPin.matches("[0-9a-f]{20}")){throw new IllegalArgumentException("Invalid member pin: "+id);}
				final long expectedResidues=tails ? lp.parseLong(residuesColumn) : -1;
				final FamilyResult result=buildFamily(fasta,geometry,seed,widths,registerBits,exponent,positionMode,partitions,overlap,positions);
				if(result.members()!=expected || (tails && result.residues!=expectedResidues)){
					throw new IllegalStateException("Member/residue count mismatch: "+id);
				}
				requirePin(fasta,fastaPin,hashBuffer);
				for(int w=0; w<widths.length; w++){
					for(boolean includeEmpty : new boolean[]{true,false}){
						lineOut.clear();
						lineOut.append(namespace).tab().append(id).tab().append(widths[w]).tab()
							.append(expected).tab().append(result.residues).tab().append(result.words).tab()
							.append(includeEmpty ? "literal" : "nonempty").tab();
						for(int kind=0; kind<3; kind++){
							if(kind>0){lineOut.tab();}
							if(registerBits==16){appendArray(lineOut,result.profiles[w],includeEmpty,kind);}
							else{appendArray(lineOut,result.wideProfiles[w],includeEmpty,kind);}
						}
						if(result.positions!=null){
							lineOut.tab();
							appendPositions(lineOut,result.profiles[w],result.positions[w],includeEmpty);
						}
						if(registerBits==32){lineOut.append("\t2\t32\t").append(k).tab().append(seed).tab().append(exponent).tab().append(ENCODING32);}
						lineOut.nl();
						profiles.print(lineOut);
					}
				}
				families++; members+=expected; residues+=result.residues; words+=result.words;
				System.err.println("family_complete\t"+namespace+"\t"+id+"\t"+expected+"\t"+families);
			}
			complete=true;
		}finally{
			final boolean inputError=table.close(), outputError=profiles.poisonAndWait();
			if(inputError || outputError){throw new IllegalStateException("Manifest/profile I/O failed; no PASS");}
		}
		if(!complete || families==0){throw new IllegalStateException("No families built");}
		requirePin(manifestPath,pin,hashBuffer);
		final boolean legacyK5=geometry.isExact(5);
		final String summary="namespace\tmanifest_sha80\tk\tseed\texponent\tshard\tshards\tmanifest_families\tbuilt_families\tmembers\tresidues\tvalid_words\telapsed_seconds\tsketch"+
			(legacyK5 ? "" : "\tgeometry")+
			(registerBits==32 ? "\tversion\tregisterbits\tencoding" : "")+
			positionMetadataColumns(positionMode)+(official==null ? "" : "\tblacklist_sha80\tblacklist_words\t"+official.metadataColumns())+"\n"+
			namespace+"\t"+pin+"\t"+k+"\t"+seed+"\t"+exponent+"\t"+shard+"\t"+shards+"\t"+row+"\t"+families+"\t"+members+"\t"+residues+"\t"+words+"\t"+((System.nanoTime()-started)*1e-9)+
			"\t"+sketch+(legacyK5 ? "" : "\t"+geometry.name())+(registerBits==32 ? "\t2\t32\t"+ENCODING32 : "")+
			positionMetadataValues(positionMode,partitions,overlap,geometry)+(official==null ? "" : "\t"+official.sha80()+"\t"+official.size()+"\t"+official.metadataValues())+"\n";
		Files.write(output.resolve("summary.tsv"),summary.getBytes(StandardCharsets.UTF_8));
		Files.write(output.resolve("recipe.txt"),String.join("\n",args).concat("\n").getBytes(StandardCharsets.UTF_8));
		Files.write(output.resolve("PASS"),"DDL_FAMILY_PROFILE_BUILD_PASS\n".getBytes(StandardCharsets.US_ASCII));
	}

	private static int column(String[] header,String name){
		final int index=Arrays.asList(header).indexOf(name);
		if(index<0){throw new IllegalArgumentException("Missing manifest column "+name);}
		return index;
	}

	/** Ascending unique powers of two keep per-member folding unambiguous. */
	private static int[] parseWidths(String value){
		final String[] parts=value.split(",",-1);
		final int[] widths=new int[parts.length];
		int previous=0;
		for(int i=0; i<parts.length; i++){
			final int width=Integer.parseInt(parts[i]);
			if(width<128 || width>MAX_BUCKETS || Integer.bitCount(width)!=1 || width<=previous){
				throw new IllegalArgumentException("widths must be ascending unique powers of two,128.."+MAX_BUCKETS);
			}
			widths[i]=previous=width;
		}
		assert(widths.length>0) : "A build needs at least one profile width";
		return widths;
	}

	private static FamilyResult buildFamily(Path path,ProteinKmerGeometry geometry,long seed,int[] widths,int registerBits,int exponent,ProteinPositionMode positionMode,int partitions,int overlap,boolean positions){
		assert(widths.length>0) : "The largest requested member sketch supplies all smaller widths";
		final FamilyResult result=new FamilyResult(widths,registerBits);
		final char[][] sketches=registerBits==16 ? new char[widths.length][] : null;
		final int[][] wideSketches=registerBits==32 ? new int[widths.length][] : null;
		for(int i=0; i<widths.length; i++){
			if(registerBits==16){sketches[i]=new char[widths[i]];}else{wideSketches[i]=new int[widths[i]];}
		}
		final int largest=widths.length-1;
		final FileFormat format=FileFormat.testInput(path.toString(),FileFormat.FASTA,null,false,true);
		final Streamer reader=StreamerFactory.makeStreamer(format,null,true,-1,false,true,1);
		reader.start();
		try{
			for(ListNum<Read> batch=reader.nextList(); batch!=null; batch=reader.nextList()){
				for(Read read : batch.list){
					if(read==null || read.bases==null || read.mate!=null){throw new IllegalArgumentException("Require unpaired protein members");}
					result.residues+=read.bases.length;
					if(registerBits==16){
						result.words+=positionMode.positional() ? sketchPositional(read.bases,geometry,seed,partitions,overlap,positionMode,sketches[largest],null) :
							sketch(read.bases,geometry,seed,sketches[largest],null);
						for(int i=0; i<largest; i++){fold(sketches[largest],sketches[i]);}
						for(int i=0; i<widths.length; i++){result.profiles[i].add(sketches[i]);}
					}else{
						result.words+=sketch(read.bases,geometry.exact(),seed,exponent,wideSketches[largest]);
						for(int i=0; i<largest; i++){fold(wideSketches[largest],wideSketches[i]);}
						for(int i=0; i<widths.length; i++){result.wideProfiles[i].add(wideSketches[i]);}
					}
				}
			}
		}finally{
			if(ReadWrite.closeStream(reader)){throw new IllegalStateException("Failed reading protein members: "+path);}
		}
		assert(result.members()==(registerBits==16 ? result.profiles[largest].members() : result.wideProfiles[largest].members())) :
			"Every member must vote once at every width";
		if(positions){addPositions(path,geometry,seed,partitions,overlap,result,positionMode);}
		return result;
	}

	/** Exact second pass: only modal-register members vote for stop-relative coordinates. */
	private static void addPositions(Path path,ProteinKmerGeometry geometry,long seed,int partitions,int overlap,FamilyResult result,ProteinPositionMode positionMode){
		assert(result.profiles!=null && result.profiles.length==1 && result.profiles[0].buckets()==POSITIONAL_BUCKETS) :
			"Position payloads are currently bound to one positional 1024-bucket profile";
		final char[] sketch=new char[POSITIONAL_BUCKETS];
		final DDLQueryOrigins origins=new DDLQueryOrigins(POSITIONAL_BUCKETS);
		final PositionProfile profile=new PositionProfile(POSITIONAL_BUCKETS);
		final FileFormat format=FileFormat.testInput(path.toString(),FileFormat.FASTA,null,false,true);
		final Streamer reader=StreamerFactory.makeStreamer(format,null,true,-1,false,true,1);
		reader.start();
		int members=0;
		try{
			for(ListNum<Read> batch=reader.nextList(); batch!=null; batch=reader.nextList()){
				for(Read read : batch.list){
					if(read==null || read.bases==null || read.mate!=null){throw new IllegalArgumentException("Require unpaired protein members");}
					sketchPositional(read.bases,geometry,seed,partitions,overlap,sketch,null,origins,positionMode);
					profile.add(sketch,origins,result.profiles[0],read.bases.length);
					members++;
				}
			}
		}finally{
			if(ReadWrite.closeStream(reader)){throw new IllegalStateException("Failed rereading protein members: "+path);}
		}
		if(members!=result.members()){throw new IllegalStateException("Position pass member count mismatch: "+path);}
		result.positions=new PositionProfile[]{profile};
	}

	private static void appendArray(ByteBuilder out,Profile profile,boolean includeEmpty,int kind){
		assert(kind>=0 && kind<=2) : "Profile arrays are value, support, occupancy";
		for(int b=0; b<profile.buckets(); b++){
			if(b>0){out.append(',');}
			out.append(kind==0 ? (int)profile.value(b,includeEmpty) : kind==1 ? profile.votes(b,includeEmpty) : profile.occupied(b));
		}
	}

	private static void appendPositions(ByteBuilder out,Profile profile,PositionProfile positions,boolean includeEmpty){
		for(int b=0; b<profile.buckets(); b++){
			if(b>0){out.append(',');}
			out.append(positions.position(b,profile,includeEmpty));
		}
	}

	private static void appendArray(ByteBuilder out,Profile32 profile,boolean includeEmpty,int kind){
		assert(kind>=0 && kind<=2) : "Wide profile arrays preserve uint32 values and int support/occupancy";
		for(int b=0; b<profile.buckets(); b++){
			if(b>0){out.append(',');}
			out.append(kind==0 ? Integer.toUnsignedLong(profile.value(b,includeEmpty)) :
				kind==1 ? profile.votes(b,includeEmpty) : profile.occupied(b));
		}
	}

	private static void requirePin(Path path,String expected,byte[] buffer) throws Exception{
		final MessageDigest digest=MessageDigest.getInstance("SHA-256");
		try(InputStream input=Files.newInputStream(path)){
			for(int n=input.read(buffer); n!=-1; n=input.read(buffer)){digest.update(buffer,0,n);}
		}
		final byte[] bytes=digest.digest();
		final StringBuilder suffix=new StringBuilder(20);
		for(int i=bytes.length-10; i<bytes.length; i++){
			final int value=bytes[i]&255;
			suffix.append("0123456789abcdef".charAt(value>>>4)).append("0123456789abcdef".charAt(value&15));
		}
		if(!suffix.toString().equals(expected)){throw new IllegalStateException("sha80 mismatch: "+path);}
	}

	private static final class FamilyResult{
		FamilyResult(int[] widths,int registerBits){
			profiles=registerBits==16 ? new Profile[widths.length] : null;
			wideProfiles=registerBits==32 ? new Profile32[widths.length] : null;
			for(int i=0; i<widths.length; i++){
				if(registerBits==16){profiles[i]=new Profile(widths[i]);}else{wideProfiles[i]=new Profile32(widths[i]);}
			}
		}
		int members(){return profiles!=null ? profiles[0].members() : wideProfiles[0].members();}
		final Profile[] profiles;
		final Profile32[] wideProfiles;
		PositionProfile[] positions;
		long residues,words;
	}

	/** Bounds dense search offsets and keeps packed bucket/value keys in signed int. */
	public static final int MAX_BUCKETS=16384;
	public static final String ENCODING32="uint32-decimal";
	public static final String HEADER16="namespace\tfamily_id\tbuckets\tmembers\tresidues\tvalid_words\tmode\tvalues\tsupport\toccupancy";
	public static final String HEADER32=HEADER16+"\tversion\tregisterbits\tk\tseed\texponent\tencoding";
	public static final String SKETCH_PLAIN="plain";
	public static final String SKETCH_POSITIONAL="positional-v1";
	public static final int POSITIONAL_K=5;
	public static final int POSITIONAL_BUCKETS=1024;
	public static final int POSITIONAL_QUADRANTS=4;
	public static final int POSITIONAL_QUADRANT_WORDS=49;
	public static final int POSITIONAL_OVERLAP=8;
	public static final int POSITIONAL_TAIL_WORDS=196;
	public static final String HEADER16_POSITIONS=HEADER16+"\tpositions";
	public static final String RELATIVE_POSITION_BASIS="residue-length";

	/** Width changes only the bucket mask; absolute uint16 scores and positional tags stay unchanged. */
	public static void validatePositionalWidth(final int width){
		if(width!=512 && width!=POSITIONAL_BUCKETS){
			throw new IllegalArgumentException("Positional sketches require 512 or 1024 registers: "+width);
		}
	}

	/** Positional summary columns; fixed-profile headers remain byte-identical. */
	public static String positionMetadataColumns(ProteinPositionMode mode){
		if(mode==null){throw new IllegalArgumentException("Profile metadata requires a position mode");}
		return !mode.positional() ? "" : "\tpartitions\toverlap\ttag_bits\t"+
			(mode.fixedCoordinates() ? "segment_words" : "position_basis");
	}

	/** Relative layouts declare their basis rather than inventing a constant segment width. */
	public static String positionMetadataValues(ProteinPositionMode mode,int partitions,int overlap,ProteinKmerGeometry geometry){
		if(mode==null){throw new IllegalArgumentException("Profile metadata requires a position mode");}
		if(!mode.positional()){return "";}
		validatePositionalLayout(partitions,overlap,geometry,mode);
		return "\t"+partitions+"\t"+overlap+"\t"+positionalTagBits(partitions)+"\t"+
			(mode.fixedCoordinates() ? Integer.toString(positionalSegmentWords(partitions,geometry)) : RELATIVE_POSITION_BASIS);
	}

	/** The tag is part of the hashed word, so unsupported layouts must never alias an existing recipe. */
	public static int positionalTagBits(int partitions){
		if(partitions<2 || partitions>16 || Integer.bitCount(partitions)!=1){
			throw new IllegalArgumentException("Positional partitions must be 2, 4, 8 or 16: "+partitions);
		}
		return Integer.numberOfTrailingZeros(partitions);
	}

	/** Fixed stop-relative segments cover one 200-residue tail's protein-word starts, never rescaled by member length. */
	public static int positionalSegmentWords(int partitions){
		positionalTagBits(partitions);
		return positionalSegmentWords(partitions, POSITIONAL_K);
	}

	/** Stop-relative segments for one configured protein-word span within a 200-residue tail. */
	public static int positionalSegmentWords(int partitions, ProteinKmerGeometry geometry){
		if(geometry==null){throw new IllegalArgumentException("Positional geometry is required");}
		return positionalSegmentWords(partitions, geometry.span());
	}

	private static int positionalSegmentWords(int partitions,int span){
		positionalTagBits(partitions);
		final int tailWords=201-span;
		if(tailWords<1){throw new IllegalArgumentException("Positional protein word cannot exceed a 200-residue tail");}
		final int words=(tailWords+partitions-1)/partitions;
		assert(words*partitions>=tailWords && (words-1)*partitions<tailWords) :
			"Ceiling division must cover all reference tail words without adding a whole segment";
		return words;
	}

	/** Only immediate neighbours are duplicated; an overlap cannot cross a whole intervening segment. */
	public static void validatePositionalLayout(int partitions,int overlap){
		final int words=positionalSegmentWords(partitions);
		if(overlap<0 || overlap>words){throw new IllegalArgumentException("Positional overlap must be 0.."+words+": "+overlap);}
	}

	/** Validate overlap after accounting for the configured word span. */
	public static void validatePositionalLayout(int partitions,int overlap,ProteinKmerGeometry geometry){
		final int words=positionalSegmentWords(partitions,geometry);
		if(overlap<0 || overlap>words){throw new IllegalArgumentException("Positional overlap must be 0.."+words+": "+overlap);}
	}

	/** Quadrant placement reserves exactly two selector bits; other tag counts cannot share its identity. */
	public static void validatePositionalLayout(int partitions,int overlap,ProteinKmerGeometry geometry,ProteinPositionMode mode){
		if(mode==null || !mode.positional()){
			throw new IllegalArgumentException("Positional layout requires an explicit positional mode");
		}
		validatePositionalLayout(partitions,overlap,geometry);
		if(mode==ProteinPositionMode.FIXED_QUADRANT && partitions!=4){
			throw new IllegalArgumentException("Fixed-quadrant bucket placement requires exactly four partitions: "+partitions);
		}
	}

	/** Legacy encoding is fixed; wide encoding is instance-bound, never a global mutation. */
	public static void validateLayout(int registerBits,int exponent){
		if((registerBits!=16 && registerBits!=32) || exponent<1 || exponent>8 || (registerBits==16 && exponent!=6)){
			throw new IllegalArgumentException("Family profiles require registerbits=16 exponent=6 or registerbits=32 exponent=1..8");
		}
	}

	/** Reuse the production wide score encoding, without allocating a tracker per member. */
	public static long sketch(byte[] sequence,int k,long seed,int exponent,int[] out){
		validateLayout(32,exponent);
		if(sequence==null || k<1 || k>12 || out==null || out.length<1 || Integer.bitCount(out.length)!=1){
			throw new IllegalArgumentException("Invalid wide protein sketch parameters");
		}
		Arrays.fill(out,0);
		final long mask=(1L<<(5*k))-1;
		final int bucketMask=out.length-1;
		long word=0,words=0; int valid=0;
		for(byte residue : sequence){
			final int aa=residue>=0 && residue<AminoAcid.acidToNumberNoStops.length ? AminoAcid.acidToNumberNoStops[residue] : -1;
			if(aa<0){word=0; valid=0; continue;}
			word=((word<<5)|aa)&mask;
			if(++valid<k){continue;}
			final long hash=Tools.hash64shift(word^seed);
			final int bucket=(int)(hash&bucketMask),value=DynamicDemiLog32.scoreHash(hash,exponent);
			if(Integer.compareUnsigned(value,out[bucket])>0){out[bucket]=value;}
			words++;
		}
		assert(words<=Math.max(0,sequence.length-k+1)) : "Each wide word occurrence starts inside this member";
		return words;
	}

	/** Fold member maxima with unsigned ordering, before any family votes are accumulated. */
	public static void fold(int[] source,int[] target){
		if(source==null || target==null || source==target || target.length<1 ||
			Integer.bitCount(target.length)!=1 || Integer.bitCount(source.length)!=1 || target.length>source.length){
			throw new IllegalArgumentException("Invalid wide sketch fold");
		}
		Arrays.fill(target,0);
		for(int i=0,mask=target.length-1; i<source.length; i++){
			final int b=i&mask;
			if(Integer.compareUnsigned(source[i],target[b])>0){target[b]=source[i];}
		}
		assert(source.length%target.length==0) : "Complete power-of-two member bucket groups preserve direct sketch maxima";
	}

	/** Separate wide votes: a cast before shifting preserves every bucket and all value bits. */
	public static final class Profile32{
		public Profile32(int buckets){
			if(buckets<1 || buckets>MAX_BUCKETS || Integer.bitCount(buckets)!=1){throw new IllegalArgumentException("Invalid wide profile width");}
			mode=new int[buckets]; support=new int[buckets]; occupied=new int[buckets]; counts=new LongHashMap(buckets*2);
		}
		public void add(int[] member){
			if(member==null || member.length!=mode.length || members==Integer.MAX_VALUE){throw new IllegalArgumentException("Wide member width/count mismatch");}
			members++;
			for(int b=0; b<member.length; b++){
				final int value=member[b]; if(value==0){continue;}
				occupied[b]++;
				final int count=counts.increment(((long)b<<32)|Integer.toUnsignedLong(value),1);
				if(count>support[b] || (count==support[b] && Integer.compareUnsigned(value,mode[b])<0)){
					mode[b]=value; support[b]=count;
				}
				assert(count<=members) : "Each actual member votes once per bucket/full uint32 value";
			}
		}
		public int value(int bucket,boolean includeEmpty){
			assert(bucket>=0 && bucket<mode.length) : "Wide mode lookup must use a bucket in this profile";
			return includeEmpty && members-occupied[bucket]>=support[bucket] ? 0 : mode[bucket];
		}
		public int votes(int bucket,boolean includeEmpty){
			assert(bucket>=0 && bucket<mode.length) : "Wide support lookup must use a bucket in this profile";
			return includeEmpty && value(bucket,true)==0 ? members-occupied[bucket] : support[bucket];
		}
		public int members(){return members;}
		public int buckets(){return mode.length;}
		public int occupied(int bucket){return occupied[bucket];}
		private int members;
		private final int[] mode,support,occupied;
		private final LongHashMap counts;
	}

	/** Sketch forward standard-amino-acid words; ambiguity and stops break words.
	 * Returns the number of valid word occurrences, not distinct words.
	 * Caller owns and can reuse the output array across members.
	 */
	public static long sketch(byte[] sequence, int k, long seed, char[] out){
		return sketch(sequence, k, seed, out, null);
	}

	/** Exclude literal AA pentamers before selecting register winners; retain raw valid-word totals. */
	public static long sketch(byte[] sequence, int k, long seed, char[] out, DDLPentamerBlacklist blacklist){
		return sketch(sequence, ProteinKmerGeometry.exact(k), seed, out, blacklist);
	}

	/** Exclude optional exact-prefix pentamers before selecting register winners; retain raw valid-word totals. */
	public static long sketch(byte[] sequence, ProteinKmerGeometry geometry, long seed, char[] out, DDLPentamerBlacklist blacklist){
		return sketchSequence(sequence,geometry,seed,out,blacklist,false);
	}

	/** Same sketch from BBTools standard codes 0..19 and the extended X code.
	 * This prevents encoded clustering residues being mistaken for ASCII input.
	 */
	public static long sketchEncoded(byte[] sequence,int k,long seed,char[] out,DDLPentamerBlacklist blacklist){
		return sketchSequence(sequence,ProteinKmerGeometry.exact(k),seed,out,blacklist,true);
	}

	public static long sketchEncoded(byte[] sequence,ProteinKmerGeometry geometry,long seed,char[] out,DDLPentamerBlacklist blacklist){
		return sketchSequence(sequence,geometry,seed,out,blacklist,true);
	}

	/** Shares the word/hash kernel between ASCII and encoded sequence storage. */
	private static long sketchSequence(byte[] sequence,ProteinKmerGeometry geometry,long seed,char[] out,DDLPentamerBlacklist blacklist,boolean encoded){
		if(sequence==null || geometry==null || out==null || out.length<1 ||
			Integer.bitCount(out.length)!=1){throw new IllegalArgumentException("Invalid protein sketch parameters");}
		blacklist=ProteinBlacklistPolicy.resolve(geometry,blacklist);
		Arrays.fill(out, (char)0);
		final long wordMask=geometry.rollingMask();
		final int bucketMask=out.length-1;
		long word=0, validWords=0;
		int valid=0;
		for(byte residue : sequence){
			final int aa=sketchResidue(residue,encoded);
			if(aa<0){word=0; valid=0; continue;}
			word=((word<<5)|aa)&wordMask;
			if(++valid<geometry.span()){continue;}
			if(blacklist==null || !blacklist.contains(blacklist.rawWord(geometry, word))){
				store(geometry.word(word), seed, out, bucketMask);
			}
			validWords++;
		}
		assert(validWords<=Math.max(0, sequence.length-geometry.span()+1)) : "Each word must start within the member sequence";
		return validWords;
	}

	/** Legacy entry point preserves the four-partition, overlap-eight hashed words exactly. */
	public static long sketchPositional(byte[] sequence,int k,long seed,char[] out){
		return sketchPositional(sequence,k,seed,POSITIONAL_QUADRANTS,POSITIONAL_OVERLAP,out);
	}

	/** Sketch protein words with fixed stop-relative tags; short tails cannot vote into unoccupied segments. */
	public static long sketchPositional(byte[] sequence,int k,long seed,int partitions,int overlap,char[] out){
		return sketchPositional(sequence, k, seed, partitions, overlap, out, null);
	}

	/** Raw-word exclusion precedes every positional tag, including duplicated border votes. */
	public static long sketchPositional(byte[] sequence,int k,long seed,int partitions,int overlap,char[] out,DDLPentamerBlacklist blacklist){
		return sketchPositional(sequence, ProteinKmerGeometry.exact(k), seed, partitions, overlap, out, blacklist);
	}

	public static long sketchPositional(byte[] sequence,ProteinKmerGeometry geometry,long seed,int partitions,int overlap,char[] out,DDLPentamerBlacklist blacklist){
		return sketchPositional(sequence,geometry,seed,partitions,overlap,ProteinPositionMode.FIXED,out,blacklist);
	}

	/** Explicit fixed or relative coordinates; existing overloads retain fixed geometry. */
	public static long sketchPositional(byte[] sequence,ProteinKmerGeometry geometry,long seed,int partitions,int overlap,
			ProteinPositionMode mode,char[] out,DDLPentamerBlacklist blacklist){
		return sketchPositionalSequence(sequence,geometry,seed,partitions,overlap,out,blacklist,false,mode);
	}

	/** Encoded-residue counterpart of the exact same stop-relative positional kernel. */
	public static long sketchPositionalEncoded(byte[] sequence,int k,long seed,int partitions,int overlap,char[] out,DDLPentamerBlacklist blacklist){
		return sketchPositionalEncoded(sequence,ProteinKmerGeometry.exact(k),seed,partitions,overlap,out,blacklist);
	}

	public static long sketchPositionalEncoded(byte[] sequence,ProteinKmerGeometry geometry,long seed,int partitions,int overlap,char[] out,DDLPentamerBlacklist blacklist){
		return sketchPositionalEncoded(sequence,geometry,seed,partitions,overlap,out,blacklist,ProteinPositionMode.FIXED);
	}

	/** Encoded members and raw queries share the same explicit coordinate mode. */
	public static long sketchPositionalEncoded(byte[] sequence,ProteinKmerGeometry geometry,long seed,int partitions,int overlap,char[] out,
			DDLPentamerBlacklist blacklist,ProteinPositionMode mode){
		return sketchPositionalSequence(sequence,geometry,seed,partitions,overlap,out,blacklist,true,mode);
	}

	/** Shared positional loop; relative mode normalizes stop-relative distances by residue length. */
	private static long sketchPositionalSequence(byte[] sequence,ProteinKmerGeometry geometry,long seed,int partitions,int overlap,char[] out,
			DDLPentamerBlacklist blacklist,boolean encoded,ProteinPositionMode mode){
		if(sequence==null || geometry==null || out==null || mode==null || !mode.positional()){
			throw new IllegalArgumentException("Positional sketches require sequence, geometry, output and an explicit fixed or relative mode");
		}
		validatePositionalWidth(out.length);
		validatePositionalLayout(partitions,overlap,geometry,mode);
		blacklist=ProteinBlacklistPolicy.resolve(geometry,blacklist);
		Arrays.fill(out, (char)0);
		final int lastStart=sequence.length-geometry.span();
		if(lastStart<0){return 0;}
		final int segmentWords=positionalSegmentWords(partitions,geometry), tagBits=positionalTagBits(partitions);
		final int maxQuadrant=positionTag(lastStart,sequence.length,partitions,segmentWords,mode);
		final long wordMask=geometry.rollingMask();
		long word=0,validWords=0;
		int valid=0;
		for(int pos=0; pos<sequence.length; pos++){
			final byte residue=sequence[pos];
			final int aa=sketchResidue(residue,encoded);
			if(aa<0){word=0; valid=0; continue;}
			word=((word<<5)|aa)&wordMask;
			if(++valid<geometry.span()){continue;}
			final int ordinal=lastStart-(pos-geometry.span()+1);
			if(blacklist==null || !blacklist.contains(blacklist.rawWord(geometry, word))){
				storePositional(geometry.word(word),seed,ordinal,maxQuadrant,partitions,segmentWords,tagBits,overlap,out,out.length-1,
					sequence.length,mode);
			}
			validWords++;
		}
		assert(validWords<=Math.max(0, sequence.length-geometry.span()+1)) : "Each positional word must start within the member sequence";
		return validWords;
	}

	/** Enumerate every filtered encoded key, before per-bucket maximum reduction.
	 * Shares residue decoding, geometry and position boundaries with ordinary sketches.
	 * This is consensus membership for tie-breaking, not extra family support.
	 */
	public static long enumerateKeys(byte[] sequence,ProteinKmerGeometry geometry,long seed,
			ProteinPositionMode mode,int partitions,int overlap,DDLPentamerBlacklist blacklist,
			DDLKeyTieBreak.KeySet output){
		return enumerateKeys(sequence,geometry,seed,mode,partitions,overlap,blacklist,output,false);
	}

	/** Encoded counterpart; filtering, boundaries and tags share the ASCII kernel. */
	public static long enumerateKeysEncoded(byte[] sequence,ProteinKmerGeometry geometry,long seed,
			ProteinPositionMode mode,int partitions,int overlap,DDLPentamerBlacklist blacklist,
			DDLKeyTieBreak.KeySet output){
		return enumerateKeys(sequence,geometry,seed,mode,partitions,overlap,blacklist,output,true);
	}

	private static long enumerateKeys(byte[] sequence,ProteinKmerGeometry geometry,long seed,
			ProteinPositionMode mode,int partitions,int overlap,DDLPentamerBlacklist blacklist,
			DDLKeyTieBreak.KeySet output,boolean encoded){
		if(sequence==null || geometry==null || mode==null || output==null){
			throw new IllegalArgumentException("All-key enumeration requires sequence, geometry, mode and output");
		}
		if(mode.positional()){validatePositionalWidth(output.width());validatePositionalLayout(partitions,overlap,geometry,mode);}
		blacklist=ProteinBlacklistPolicy.resolve(geometry,blacklist);output.clear();
		final int lastStart=sequence.length-geometry.span();
		if(lastStart<0){return 0;}
		final int segment=mode.positional() ? positionalSegmentWords(partitions,geometry) : 0;
		final int tagBits=mode.positional() ? positionalTagBits(partitions) : 0;
		final int maxTag=mode.positional() ? positionTag(lastStart,sequence.length,partitions,segment,mode) : 0;
		final long mask=geometry.rollingMask();long word=0,validWords=0;int valid=0;
		for(int pos=0;pos<sequence.length;pos++){
			final int aa=sketchResidue(sequence[pos],encoded);
			if(aa<0){word=0;valid=0;continue;}
			word=((word<<5)|aa)&mask;
			if(++valid<geometry.span()){continue;}
			validWords++;
			if(blacklist!=null && blacklist.contains(blacklist.rawWord(geometry,word))){continue;}
			final long encodedWord=geometry.word(word);
			if(!mode.positional()){observeKey(encodedWord,seed,output,0,mode);continue;}
			final int ordinal=lastStart-(pos-geometry.span()+1);
			final int tag=positionTag(ordinal,sequence.length,partitions,segment,mode);
			observeKey((encodedWord<<tagBits)|tag,seed,output,tag,mode);
			if(ordinal-positionBoundary(tag,sequence.length,partitions,segment,mode)<overlap && tag>0){
				observeKey((encodedWord<<tagBits)|(tag-1),seed,output,tag-1,mode);
			}
			if(positionBoundary(tag+1,sequence.length,partitions,segment,mode)-ordinal<=overlap && tag<maxTag){
				observeKey((encodedWord<<tagBits)|(tag+1),seed,output,tag+1,mode);
			}
		}
		return validWords;
	}

	private static void observeKey(long word,long seed,DDLKeyTieBreak.KeySet output,int tag,ProteinPositionMode mode){
		final long hash=Tools.hash64shift(word^seed);
		final int nlz=Long.numberOfLeadingZeros(hash)&63;
		final int score=Math.max(1,(nlz<<10)+(int)((~(hash>>>(53-nlz)))&1023));
		output.add(sketchBucket(hash,output.width()-1,tag,mode),score);
	}

	/** Raw ambiguity breaks words; encoded input follows Blosum62.encode (0..19 or X).
	 * A stop or ASCII byte in encoded storage is malformed, not an empty sketch.
	 */
	private static int sketchResidue(byte residue,boolean encoded){
		if(!encoded){
			return residue>=0 && residue<AminoAcid.acidToNumberNoStops.length ? AminoAcid.acidToNumberNoStops[residue] : -1;
		}
		if(residue>=0 && residue<20){return residue;}
		if(residue==AminoAcid.acidToNumberExtended['X']){return -1;}
		throw new IllegalArgumentException("Expected encoded standard amino acid or X, got "+residue);
	}

	/** Diagnostic variant retaining residue starts of each final register winner.
	 * The ordinary kernel above remains untouched; parity fixtures compare its bytes.
	 */
	public static long sketchPositional(byte[] sequence,int k,long seed,int partitions,int overlap,char[] out,
			DDLPentamerBlacklist blacklist,DDLQueryOrigins origins){
		return sketchPositional(sequence, ProteinKmerGeometry.exact(k), seed, partitions, overlap, out, blacklist, origins);
	}

	public static long sketchPositional(byte[] sequence,ProteinKmerGeometry geometry,long seed,int partitions,int overlap,char[] out,
			DDLPentamerBlacklist blacklist,DDLQueryOrigins origins){
		return sketchPositional(sequence,geometry,seed,partitions,overlap,out,blacklist,origins,ProteinPositionMode.FIXED);
	}

	/** Origin tracing uses the same primary tags and integer overlap boundaries as ordinary sketches. */
	public static long sketchPositional(byte[] sequence,ProteinKmerGeometry geometry,long seed,int partitions,int overlap,char[] out,
			DDLPentamerBlacklist blacklist,DDLQueryOrigins origins,ProteinPositionMode mode){
		if(sequence==null || geometry==null || out==null || origins==null || mode==null || !mode.positional()){
			throw new IllegalArgumentException("Origin tracing requires positional sketches and reusable origin scratch");
		}
		validatePositionalWidth(out.length);
		validatePositionalLayout(partitions,overlap,geometry,mode);
		blacklist=ProteinBlacklistPolicy.resolve(geometry,blacklist);
		Arrays.fill(out,(char)0); origins.reset(out.length);
		final int lastStart=sequence.length-geometry.span();
		if(lastStart<0){return 0;}
		final int segmentWords=positionalSegmentWords(partitions,geometry),tagBits=positionalTagBits(partitions);
		final int maxQuadrant=positionTag(lastStart,sequence.length,partitions,segmentWords,mode);
		final long wordMask=geometry.rollingMask();
		long word=0,validWords=0;
		int valid=0;
		for(int pos=0; pos<sequence.length; pos++){
			final byte residue=sequence[pos];
			final int aa=residue>=0 && residue<AminoAcid.acidToNumberNoStops.length ?
				AminoAcid.acidToNumberNoStops[residue] : -1;
			if(aa<0){word=0; valid=0; continue;}
			word=((word<<5)|aa)&wordMask;
			if(++valid<geometry.span()){continue;}
			final int start=pos-geometry.span()+1,ordinal=lastStart-start;
			if(blacklist==null || !blacklist.contains(blacklist.rawWord(geometry, word))){
				final int q=positionTag(ordinal,sequence.length,partitions,segmentWords,mode);
				final long lower=positionBoundary(q,sequence.length,partitions,segmentWords,mode);
				final long upper=positionBoundary(q+1,sequence.length,partitions,segmentWords,mode);
				assert(q<=maxQuadrant) : "Tail length clips the occupied stop-relative tags";
				final long hashWord=geometry.word(word);
				storeOrigin((hashWord<<tagBits)|q,seed,out,origins,start,q,mode);
				if(ordinal-lower<overlap && q>0){storeOrigin((hashWord<<tagBits)|(q-1),seed,out,origins,start,q-1,mode);}
				if(upper-ordinal<=overlap && q<maxQuadrant){storeOrigin((hashWord<<tagBits)|(q+1),seed,out,origins,start,q+1,mode);}
			}
			validWords++;
		}
		assert(validWords<=Math.max(0,sequence.length-geometry.span()+1)) : "Each tagged vote has an in-sequence protein word start";
		return validWords;
	}

	/** Same absolute DDL encoding as store(), recording origins before replacing the winner. */
	private static void storeOrigin(long word,long seed,char[] out,DDLQueryOrigins origins,int start,int tag,ProteinPositionMode mode){
		final long key=Tools.hash64shift(word^seed);
		final int nlz=Long.numberOfLeadingZeros(key)&63;
		final int score=Math.max(1,(nlz<<10)+(int)((~(key>>>(53-nlz)))&1023));
		final int bucket=sketchBucket(key,out.length-1,tag,mode);
		assert(score>0 && score<=65535) : "Origin-traced DDL values use the same unsigned16 encoding as ordinary sketching";
		origins.observe(bucket,score,out[bucket],start,tag);
		if(score>out[bucket]){out[bucket]=(char)score;}
	}

	private static void storePositional(long word,long seed,int ordinal,int maxQuadrant,int partitions,
			int segmentWords,int tagBits,int overlap,char[] out,int bucketMask,int length,ProteinPositionMode mode){
		assert(ordinal>=0) : "Stop-relative ordinals are derived only from valid in-sequence protein-word starts";
		assert(maxQuadrant>=0 && maxQuadrant<partitions) : "Tail length clips the highest populated segment";
		final int q=positionTag(ordinal,length,partitions,segmentWords,mode);
		assert(q<=maxQuadrant) : "A valid word ordinal cannot exceed the clipped maximum quadrant";
		store((word<<tagBits)|q,seed,out,bucketMask,q,mode);
		final long lower=positionBoundary(q,length,partitions,segmentWords,mode);
		final long upper=positionBoundary(q+1,length,partitions,segmentWords,mode);
		if(ordinal-lower<overlap && q>0){store((word<<tagBits)|(q-1),seed,out,bucketMask,q-1,mode);}
		if(upper-ordinal<=overlap && q<maxQuadrant){store((word<<tagBits)|(q+1),seed,out,bucketMask,q+1,mode);}
	}

	/** Residue-coordinate primary tag. Long arithmetic preserves long-protein distances. */
	private static int positionTag(int ordinal,int length,int partitions,int segmentWords,ProteinPositionMode mode){
		assert(ordinal>=0 && ordinal<length && mode!=null && mode.positional()) :
			"Positional tags require an in-sequence stop distance and a positional mode";
		return mode==ProteinPositionMode.RELATIVE ? (int)((long)ordinal*partitions/length) :
			Tools.min(partitions-1,ordinal/segmentWords);
	}

	/** First integer distance in a band: ceiling(boundary*length/partitions) for relative mode.
	 * Overlap still duplicates only immediate occupied neighbors, in residue-distance units.
	 */
	private static long positionBoundary(int boundary,int length,int partitions,int segmentWords,ProteinPositionMode mode){
		assert(boundary>=0 && boundary<=partitions) : "A positional border must bound one configured band";
		return mode==ProteinPositionMode.RELATIVE ? ((long)boundary*length+partitions-1)/partitions : (long)boundary*segmentWords;
	}

	/** Same absolute 6e/10m encoding as DynamicDemiLog; no cardinality machinery.
	 * Reusing arrays avoids allocating a DDL object for every short member.
	 */
	private static void store(long word, long seed, char[] out, int bucketMask){
		store(word,seed,out,bucketMask,0,ProteinPositionMode.PLAIN);
	}

	/** Hash the original tagged word; placement never changes the hash used for the matching key. */
	private static void store(long word,long seed,char[] out,int bucketMask,int tag,ProteinPositionMode mode){
		final long key=Tools.hash64shift(word^seed);
		final int nlz=Long.numberOfLeadingZeros(key)&63;
		final int score=Math.max(1, (nlz<<10)+(int)((~(key>>>(53-nlz)))&1023));
		final int bucket=sketchBucket(key,bucketMask,tag,mode);
		assert(score>0 && score<=65535) : "Absolute DDL score must fit one unsigned 16-bit register";
		if(score>out[bucket]){out[bucket]=(char)score;}
	}

	/** Shared by winners, origins and all-key enumeration: replace only the top two selector bits. */
	private static int sketchBucket(long hash,int bucketMask,int tag,ProteinPositionMode mode){
		if(mode!=ProteinPositionMode.FIXED_QUADRANT){return (int)(hash&bucketMask);}
		assert(tag>=0 && tag<4 && (bucketMask==511 || bucketMask==1023)) :
			"Fixed-quadrant placement requires four tags and a validated 512/1024-register width; tag="+tag+" mask="+bucketMask;
		return ((int)hash&(bucketMask>>>2)) | tag*((bucketMask+1)>>>2);
	}

	/** Fold a member sketch before voting, NEVER fold a modal family profile. */
	public static void fold(char[] source, char[] target){
		if(source==null || target==null || source==target || target.length<1 ||
			Integer.bitCount(target.length)!=1 || Integer.bitCount(source.length)!=1 ||
			target.length>source.length){throw new IllegalArgumentException("Invalid sketch fold");}
		Arrays.fill(target, (char)0);
		final int mask=target.length-1;
		for(int i=0; i<source.length; i++){
			final int bucket=i&mask;
			if(source[i]>target[bucket]){target[bucket]=source[i];}
		}
		assert(source.length%target.length==0) : "Power-of-two member sketches must fold into complete bucket groups";
	}

	/** One vote per member per slot, independent of within-member word repetition.
	 * Counts use a primitive map keyed by (bucket<<16)|value. Empty votes are
	 * inferred from member count and occupancy, not stored in the hash map.
	 */
	public static final class Profile{
		public Profile(int buckets){
			this(buckets,null);
		}

		/** Optional frozen baseline depths affect only equal-support nonzero keys.
		 * Null preserves the historical smallest-unsigned-value tie rule.
		 */
		public Profile(int buckets, DDLKeyRarity rarity_){
			if(buckets<1 || buckets>MAX_BUCKETS || Integer.bitCount(buckets)!=1){
				throw new IllegalArgumentException("Profile buckets must be a power of two up to"+MAX_BUCKETS);
			}
			if(rarity_!=null && rarity_.width()!=buckets){
				throw new IllegalArgumentException("Rarity baseline must use the same profile width");
			}
			rarity=rarity_;
			mode=new char[buckets];
			support=new int[buckets];
			occupied=new int[buckets];
			counts=new IntHashMap(buckets*2);
		}

		/** Reuses one modal accumulator across family rebuilds without retaining old votes. */
		public void clear(){
			members=0;
			Arrays.fill(mode,(char)0); Arrays.fill(support,0); Arrays.fill(occupied,0);
			counts.clear();
			assert(counts.isEmpty()) : "A fresh family profile must not inherit another family bucket/value votes";
		}

		public void add(char[] member){
			if(member==null || member.length!=mode.length || members==Integer.MAX_VALUE){
				throw new IllegalArgumentException("Member width/count cannot be represented by this profile");
			}
			members++;
			for(int b=0; b<member.length; b++){
				final char value=member[b];
				if(value==0){continue;}
				occupied[b]++;
				final int count=counts.increment((b<<16)|value);
				if(count>support[b] || (count==support[b] &&
					(rarity==null ? value<mode[b] : rarity.prefers(b,value,mode[b])))){
					mode[b]=value;
					support[b]=count;
				}
				assert(count<=members) : "One member contributes at most one vote to each bucket/value";
			}
		}

		/** Nonzero ties follow the configured rule; zero retains precedence when included. */
		public char value(int bucket, boolean includeEmpty){
			assert(bucket>=0 && bucket<mode.length) : "Mode lookup requires a bucket in this profile";
			return includeEmpty && members-occupied[bucket]>=support[bucket] ? 0 : mode[bucket];
		}

		public int votes(int bucket, boolean includeEmpty){
			assert(bucket>=0 && bucket<mode.length) : "Vote lookup requires a bucket in this profile";
			return includeEmpty && value(bucket, true)==0 ? members-occupied[bucket] : support[bucket];
		}

		/** Opt-in packed selection; leaves the legacy streaming mode/support untouched.
		 * One stable pivot sort per bucket/policy, with canonical ascending input keys.
		 */
		public char[] packedValues(DDLKeyTieBreak.Policy policy,DDLKeyTieBreak.RawTailCounts raw,
				long familySeed,DDLKeyTieBreak.KeySet consensus,DDLKeyTieBreak.Workspace scratch){
			return packedValues(policy,raw,familySeed,consensus,scratch,null);
		}

		/** Optional read-only diagnostics do not participate in selection. */
		public char[] packedValues(DDLKeyTieBreak.Policy policy,DDLKeyTieBreak.RawTailCounts raw,
				long familySeed,DDLKeyTieBreak.KeySet consensus,DDLKeyTieBreak.Workspace scratch,DDLKeyTieBreak.Diagnostics diagnostics){
			if(rarity!=null || policy==null || raw==null || raw.width()!=mode.length || scratch==null ||
				scratch.candidates.length!=mode.length || consensus==null || consensus.width()!=mode.length){
				throw new IllegalArgumentException("Packed selection requires matching frozen counts, consensus and workspace");
			}
			scratch.collect(counts);final char[] result=new char[mode.length];
			for(int bucket=0;bucket<mode.length;bucket++){
				final java.util.BitSet candidates=scratch.candidates[bucket];
				final int size=candidates.cardinality();scratch.capacity(size);int i=0;
				for(int key=candidates.nextSetBit(1);key>=0;key=candidates.nextSetBit(key+1)){
					final int votes=counts.get((bucket<<16)|key);
					assert(votes>0 && votes<=members) : "Every candidate has actual member support";
					final int signal=DDLKeyTieBreak.tieSignal(policy,key,raw.count(bucket,key),familySeed,consensus.contains(bucket,key));
					scratch.keys[i]=(char)key;scratch.pivots[i]=DDLKeyTieBreak.pack(votes,signal);i++;
				}
				assert(i==size) : "No zero key in member candidate map";
				DDLKeyTieBreak.sortDescending(scratch.pivots,scratch.keys,size,scratch.pivotScratch,scratch.keyScratch,scratch.bins);
				if(size>0){
					result[bucket]=scratch.keys[0];
					assert(DDLKeyTieBreak.votes(scratch.pivots[0])==support[bucket]) : "Tie policy cannot change maximal support";
				}
				assert(policy!=DDLKeyTieBreak.Policy.SMALLEST || result[bucket]==mode[bucket]) : "Packed default preserves legacy mode";
				if(diagnostics!=null){diagnostics.observe(scratch.pivots,scratch.keys,size);}
			}
			return result;
		}

		/** Actual nonempty-key member support, distinct from any transient bonus-adjusted selection score. */
		public int memberSupport(int bucket,int key){
			if(bucket<0 || bucket>=mode.length || key<1 || key>65535){
				throw new IllegalArgumentException("Member-support lookup requires a real nonzero key in this profile");
			}
			return Math.max(0,counts.get((bucket<<16)|key));
		}

		/** Live selection into caller-owned output; immutable context, worker-owned scratch.
		 * Positive W applies once to ALL distinct consensus keys in each bin, including
		 * consensus-only keys with zero member support. Real counts/support/occupancy remain unchanged.
		 */
		public void writePackedValues(DDLKeyTieBreak.PolicyContext context,long familySeed,
				DDLKeyTieBreak.KeySet consensus,DDLKeyTieBreak.Workspace scratch,
				DDLKeyTieBreak.Diagnostics diagnostics,char[] out){
			if(rarity!=null || context==null || context.width!=mode.length || scratch==null || scratch.candidates.length!=mode.length ||
				out==null || out.length!=mode.length ||
				context.policy==DDLKeyTieBreak.Policy.RAW_TAIL_RARITY && context.rawCounts.width()!=mode.length ||
				(context.policy==DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY || context.consensusMemberBonus>0) &&
				(consensus==null || consensus.width()!=mode.length)){
				throw new IllegalArgumentException("Packed live selection requires matching policy and worker scratch");
			}
			if(context.consensusMemberBonus>0){scratch.collect(counts,consensus);}
			else{scratch.collect(counts);}
			Arrays.fill(out,(char)0);
			for(int bucket=0;bucket<mode.length;bucket++){
				final java.util.BitSet candidates=scratch.candidates[bucket];
				final int size=candidates.cardinality();scratch.capacity(size);int i=0;
				for(int key=candidates.nextSetBit(1);key>=0;key=candidates.nextSetBit(key+1)){
					final int votes=Math.max(0,counts.get((bucket<<16)|key));
					final boolean inConsensus=consensus!=null && consensus.contains(bucket,key);
					assert(votes<=members && (votes>0 || context.consensusMemberBonus>0 && inConsensus)) :
						"Only a positive bonus admits consensus-only keys; member votes remain actual observations";
					if(context.policy==DDLKeyTieBreak.Policy.RAW_TAIL_RARITY && context.rawCounts.count(bucket,key)<votes){
						throw new IllegalArgumentException("Incomplete raw view: bucket="+bucket+" key="+key+
							" member_votes="+votes+" raw_count="+context.rawCounts.count(bucket,key));
					}
					final int signal=context.policy==DDLKeyTieBreak.Policy.RAW_TAIL_RARITY ? context.rawCounts.raritySignal(bucket,key) :
						DDLKeyTieBreak.tieSignal(context.policy,key,0,familySeed,
							context.policy==DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY && inConsensus);
					final long selectionScore=Math.addExact((long)votes,inConsensus ? context.consensusMemberBonus : 0L);
					scratch.keys[i]=(char)key;scratch.pivots[i]=DDLKeyTieBreak.pack(selectionScore,signal);i++;
				}
				DDLKeyTieBreak.sortDescending(scratch.pivots,scratch.keys,size,scratch.pivotScratch,scratch.keyScratch,scratch.bins);
				if(size>0){out[bucket]=scratch.keys[0];}
				assert(context.consensusMemberBonus!=0 || context.policy!=DDLKeyTieBreak.Policy.SMALLEST || out[bucket]==mode[bucket]) :
					"W0 default mode must preserve the legacy streaming winner";
				if(diagnostics!=null){diagnostics.observe(scratch.pivots,scratch.keys,size);}
			}
		}

		public int members(){return members;}
		public int buckets(){return mode.length;}
		public int occupied(int bucket){return occupied[bucket];}
		private int members;
		private final char[] mode;
		private final int[] support, occupied;
		private final IntHashMap counts;
		private final DDLKeyRarity rarity;
	}

	/** Modal 1-based stop-relative protein-word start for each final nonzero modal register. */
	public static final class PositionProfile{
		public PositionProfile(int buckets){
			if(buckets<1 || buckets>MAX_BUCKETS || Integer.bitCount(buckets)!=1){
				throw new IllegalArgumentException("Position buckets must be a power of two up to "+MAX_BUCKETS);
			}
			positions=new int[buckets];
			support=new int[buckets];
			counts=new IntHashMap(buckets*2);
		}

		public void add(char[] member,DDLQueryOrigins origins,Profile profile,int sequenceLength){
			//TODO: Probable bug - addPositions accepts geometry.span()!=5, but this coordinate conversion assumes K5.
			// The released 512-bin study uses exact K5; repair other spans with their own position-payload fixtures.
			if(member==null || origins==null || profile==null || member.length!=positions.length || origins.first.length!=positions.length ||
					profile.buckets()!=positions.length || sequenceLength<POSITIONAL_K){
				throw new IllegalArgumentException("Position/modal profile width mismatch");
			}
			final int lastStart=sequenceLength-POSITIONAL_K;
			for(int b=0; b<member.length; b++){
				final int value=member[b];
				if(value==0 || value!=profile.value(b,false)){continue;}
				final int first=origins.first[b];
				if(first<0){throw new IllegalStateException("Modal position bucket lacks a winning origin: "+b);}
				final int position=lastStart-first+1;
				if(position<1 || position>255){throw new IllegalArgumentException("Position exceeds uint8 payload at bucket "+b+": "+position);}
				final int count=counts.increment((b<<16)|position);
				if(count>support[b] || (count==support[b] && position<positions[b])){
					positions[b]=position;
					support[b]=count;
				}
				assert(count<=profile.votes(b,false)) : "Only modal-register members contribute to a position histogram";
			}
		}

		public int position(int bucket,Profile profile,boolean includeEmpty){
			assert(bucket>=0 && bucket<positions.length) : "Position lookup requires a bucket in this profile";
			return profile.value(bucket,includeEmpty)==0 ? 0 : positions[bucket];
		}

		private final int[] positions,support;
		private final IntHashMap counts;
	}
}
