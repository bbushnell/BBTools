package ddl;

import java.io.InputStream;
import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.BasicFileAttributes;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Objects;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.LineParser1;
import shared.Shared;
import simd.SIMDFamilyProfiles;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.bam.BgzfInputStream;
import structures.ByteBuilder;
import structures.ListNum;
import structures.IntList;

/** Feasibility-only exact inverted lookup of nonzero (bucket, DDL-value) pairs.
 * Scores count equal registers, NOT sequence identity or a calibrated probability.
 * Zero/zero matches never count. Worst-case common-key postings remain explicit.
 * @author Sayu
 */
public final class DDLFamilySearch{

	public static void main(String[] args) throws Exception{
		if(args.length>0 && args[0].equals("loadbench")){loadBenchmark(Arrays.copyOfRange(args,1,args.length)); return;}
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg : args){
			final String[] parts=Parse.splitOnFirst(arg,'=');
			if(parts.length!=2 || options.put(parts[0].toLowerCase(Locale.ROOT),parts[1])!=null){throw new IllegalArgumentException("Expected unique flag=value");}
		}
		final String[] allowed={"profileslist","listsha80","in","insha80","out","buckets","mode","k","aabits","aapartition","seed","top","expected","indexed","counts","maxqueries","export","writehits","passes","formats","normalize","registerbits","exponent","sketch","positional","positionmode","partitions","overlap","blacklist","blacklistsha80","blacklistk","blacklistmode","blacklistprovenance","blacklistwords"};
		for(String key : options.keySet()){
			if(!Arrays.asList(allowed).contains(key)){throw new IllegalArgumentException("Unknown flag: "+key);}
		}
		final int width=Integer.parseInt(options.containsKey("buckets") ? options.get("buckets") : "256");
		final int k=Integer.parseInt(options.containsKey("k") ? options.get("k") : "5");
		if(options.containsKey("aabits") && options.containsKey("aapartition")){
			throw new IllegalArgumentException("Use either aabits or aapartition, not both");
		}
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(k, options.containsKey("aabits") ? options.get("aabits") : options.get("aapartition"));
		prot.FamilyBlacklistBinding.applyDefaults(options);
		final int top=Integer.parseInt(options.containsKey("top") ? options.get("top") : "10");
		final long seed=Long.parseLong(options.containsKey("seed") ? options.get("seed") : "0");
		final int registerBits=Integer.parseInt(options.containsKey("registerbits") ? options.get("registerbits") : "16");
		final int exponent=Integer.parseInt(options.containsKey("exponent") ? options.get("exponent") : "6");
		// positionmode=fixed-quadrant remains an explicit Phase A lookup hook;
		// it lost -1271/-1890 ALL/40-99 raw objective and was not adopted.
		// See results/family_sketch_phase_a_20261116.
		final ProteinPositionMode positionMode=ProteinPositionMode.resolve(options.get("positionmode"),options.get("sketch"),
			options.containsKey("positional") ? Boolean.valueOf(Parse.parseBoolean(options.get("positional"))) : null,ProteinPositionMode.PLAIN);
		final String sketchMode=positionMode.sketch;
		DDLFamilyProfiles.validateLayout(registerBits,exponent);
		if(positionMode.positional() &&
			(registerBits!=16 || (width!=512 && width!=1024) || seed!=0)){
			throw new IllegalArgumentException("Positional sketch v1 search requires registerbits=16 buckets=512/1024 seed=0");
		}
		final boolean positional=positionMode.positional();
		final int partitions=Integer.parseInt(options.getOrDefault("partitions","4"));
		final int overlap=Integer.parseInt(options.getOrDefault("overlap","8"));
		if(positional){DDLFamilyProfiles.validatePositionalLayout(partitions,overlap,geometry);}
		else if(options.containsKey("partitions") || options.containsKey("overlap")){
			throw new IllegalArgumentException("Partitions/overlap require a positional sketch");
		}
		final boolean indexed=!options.containsKey("indexed") || Parse.parseBoolean(options.get("indexed"));
		final boolean counted=options.containsKey("counts") && Parse.parseBoolean(options.get("counts"));
		final boolean normalized=options.containsKey("normalize") && Parse.parseBoolean(options.get("normalize"));
		final boolean export=options.containsKey("export") && Parse.parseBoolean(options.get("export"));
		final boolean formats=options.containsKey("formats") && Parse.parseBoolean(options.get("formats"));
		if(options.containsKey("blacklist")!=options.containsKey("blacklistsha80") ||
				(options.containsKey("blacklist") && (geometry.span()<5 || registerBits!=16 || export || formats))){
			throw new IllegalArgumentException("Blacklist requires path+sha80, 5-aa-span/16-bit sketches and metadata-preserving TSV input/output");
		}
		if(!options.containsKey("blacklist") && (options.containsKey("blacklistk") ||
				options.containsKey("blacklistmode") || options.containsKey("blacklistprovenance"))){
			throw new IllegalArgumentException("Blacklist metadata flags require blacklist=");
		}
		if(!geometry.isExact() && (registerBits==32 || export || formats)){
			throw new IllegalArgumentException("Irregular amino-acid profiles currently require metadata-preserving 16-bit TSV input/output");
		}
		if(registerBits==32 && (export || formats)){throw new IllegalArgumentException("Experimental export/formats are legacy16 only; wide TSV must not be narrowed");}
		if(positional && (export || formats)){
			throw new IllegalArgumentException("Experimental export/formats are plain-sketch only until their headers carry sketch identity");
		}
		final int passes=options.containsKey("passes") ? Integer.parseInt(options.get("passes")) : 1;
		final boolean writeHits=!options.containsKey("writehits") || Parse.parseBoolean(options.get("writehits"));
		final long maxQueries=options.containsKey("maxqueries") ? Long.parseLong(options.get("maxqueries")) : -1;
		final int scale=counted && !normalized ? 65535 : 1;
		final String mode=options.containsKey("mode") ? options.get("mode") : "literal";
		if(normalized && (!counted || !mode.equals("nonempty"))){throw new IllegalArgumentException("normalize=t requires counts=t mode=nonempty");}
		if(width<128 || width>DDLFamilyProfiles.MAX_BUCKETS || Integer.bitCount(width)!=1 || k<1 || k>12 || top<1 || top>500 || passes<1 || maxQueries==0 || maxQueries< -1 ||
			(!mode.equals("literal") && !mode.equals("nonempty"))){throw new IllegalArgumentException("Invalid profile/query parameters");}
		final Path list=Paths.get(required(options,"profileslist")).toAbsolutePath();
		final Path input=Paths.get(required(options,"in")).toAbsolutePath();
		pin(list,required(options,"listsha80"));
		pin(input,required(options,"insha80"));
		final Path output=Paths.get(required(options,"out")).toAbsolutePath();
		Files.createDirectory(output);
		Shared.setThreads(1); Shared.AMINO_IN=true;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final DDLPentamerBlacklist blacklist=ProteinBlacklistPolicy.applies(geometry) ?
			prot.FamilyBlacklistBinding.load(options,prot.FamilyCalibration.VerificationMode.SHA80) : options.containsKey("blacklist") ?
			DDLPentamerBlacklist.load(Paths.get(options.get("blacklist")),required(options,"blacklistsha80"),
				Integer.parseInt(options.getOrDefault("blacklistk","5")),
				DDLPentamerBlacklist.MatchMode.parse(options.get("blacklistmode")),
				options.get("blacklistprovenance")) : null;
		final long begin=System.nanoTime();
		final Index index=registerBits==32 ? loadForSearch32(list,width,mode,k,seed,exponent,indexed,counted,normalized) :
			load(list,width,mode,geometry,seed,sketchMode,indexed,counted,normalized,export ? output.resolve("library.bin.gz") : null,formats ? output : null,false,partitions,overlap,blacklist);
		final String expected=options.containsKey("expected") ? options.get("expected") : ".";
		int expectedIndex=-1;
		if(!expected.equals(".")){
			for(int i=0; i<index.names.length; i++){if(index.names[i].equals(expected)){expectedIndex=i; break;}}
			if(expectedIndex<0){throw new IllegalArgumentException("Expected construction family absent from library");}
		}
		final double loadSeconds=(System.nanoTime()-begin)*1e-9;
		final Workspace workspace=new Workspace(index.families());
		final char[] query=registerBits==16 ? new char[width] : null;
		final int[] wideQuery=registerBits==32 ? new int[width] : null;
		final int[] best=new int[top];
		final ByteBuilder line=new ByteBuilder();
		final ByteStreamWriter writer=writeHits ? new ByteStreamWriter(output.resolve("hits.tsv.gz").toString(),false,false,false) : null;
		if(writer!=null){
			writer.start();
			writer.print("query\trank\tfamily\t"+(normalized ? "matched_support/total_support" : counted ? "frequency_score_q16" : "matching_nonzero_registers")+"\tquery_nonzero_registers\tcandidates\tpostings\n");
		}
		final FileFormat format=FileFormat.testInput(input.toString(),FileFormat.FASTA,null,false,true);
		long queries=0,residues=0,postings=0,candidates=0,noHits=0;
		long expectedTop1=0,expectedTopK=0;
		long comparedRegisters=0;
		final long[] topHistogram=new long[width+1],expectedHistogram=new long[width+1];
		final long searchStart=System.nanoTime();
		try{
		for(int pass=0; pass<passes; pass++){
			final Streamer reader=StreamerFactory.makeStreamer(format,null,true,maxQueries,false,true,1);
			reader.start();
			try{
			for(ListNum<Read> batch=reader.nextList(); batch!=null; batch=reader.nextList()){
				for(Read read : batch.list){
					if(read==null || read.bases==null || read.mate!=null || read.id.indexOf('\t')>=0){
						throw new IllegalArgumentException("Invalid unpaired protein query");
					}
					if(registerBits==16){
						if(positional){DDLFamilyProfiles.sketchPositional(read.bases,geometry,seed,partitions,overlap,positionMode,query,blacklist);}
						else{DDLFamilyProfiles.sketch(read.bases,geometry,seed,query,blacklist);}
						index.score(query,workspace);
					}
					else{DDLFamilyProfiles.sketch(read.bases,k,seed,exponent,wideQuery); index.score(wideQuery,workspace);}
					final int found=index.top(workspace,best);
					topHistogram[found==0 ? 0 : index.histogramBin(workspace,best[0],scale)]++;
					if(expectedIndex>=0){
						expectedHistogram[index.histogramBin(workspace,expectedIndex,scale)]++;
						if(found>0 && best[0]==expectedIndex){expectedTop1++;}
						for(int rank=0; rank<found; rank++){if(best[rank]==expectedIndex){expectedTopK++;}}
					}
					int occupied=0;
					if(registerBits==16){for(char value : query){if(value>0){occupied++;}}}
					else{for(int value : wideQuery){if(value!=0){occupied++;}}}
					for(int rank=0; writer!=null && rank<Math.max(1,found); rank++){
						line.clear(); line.append(read.id).tab().append(found==0 ? 0 : rank+1).tab()
							.append(found==0 ? "." : index.names[best[rank]]).tab();
						if(normalized){line.append(found==0 ? 0 : workspace.rawScores[best[rank]]).append('/').append(found==0 ? 1 : index.totalSupport[best[rank]]);}
						else{line.append(found==0 ? 0 : workspace.scores[best[rank]]);}
						line.tab().append(occupied).tab()
							.append(workspace.size).tab().append(workspace.postings).nl();
						writer.print(line);
					}
					queries++; residues+=read.bases.length; postings+=workspace.postings; candidates+=workspace.size;
					comparedRegisters+=workspace.comparedRegisters;
					if(found==0){noHits++;}
				}
			}
			}finally{
				if(ReadWrite.closeStream(reader)){throw new IllegalStateException("Query I/O failure");}
			}
		}
		}finally{
			if(writer!=null && writer.poisonAndWait()){throw new IllegalStateException("Result I/O failure");}
		}
		final double searchSeconds=(System.nanoTime()-searchStart)*1e-9;
		pin(input,required(options,"insha80")); pin(list,required(options,"listsha80"));
		if(blacklist!=null){pin(Paths.get(options.get("blacklist")),blacklist.sha80());}
		final boolean legacyK5=geometry.isExact(5);
		final String summary="families\tbuckets\tmode\tk\tseed\tqueries\tresidues\tindex_postings\tvisited_postings\tcandidates\tno_hits\tload_seconds\tsearch_seconds\tempty_profiles\texpected_family\texpected_top1\texpected_topk\ttopk\tindexed\tcounts\tsimd\tscore_scale\tcompared_registers\tpacked_bytes\tgzip_bytes\twrite_hits\tnormalization\tsketch"+
			(legacyK5 ? "" : "\tgeometry")+
			(registerBits==32 ? "\tversion\tregisterbits\texponent\tencoding\tscan_backend" : "")+
			DDLFamilyProfiles.positionMetadataColumns(positionMode)+
			(blacklist==null ? "" : "\tblacklist_sha80\tblacklist_words\t"+blacklist.metadataColumns())+"\n"+
			index.families()+"\t"+width+"\t"+mode+"\t"+k+"\t"+seed+"\t"+queries+"\t"+residues+"\t"+
			(index.postings==null ? 0 : index.postings.length)+"\t"+postings+"\t"+candidates+"\t"+noHits+"\t"+loadSeconds+"\t"+searchSeconds+"\t"+
			index.emptyProfiles+"\t"+expected+"\t"+expectedTop1+"\t"+expectedTopK+"\t"+top+"\t"+indexed+"\t"+counted+"\t"+Shared.SIMD+"\t"+scale+"\t"+
			comparedRegisters+"\t"+index.packedBytes+"\t"+index.gzipBytes+"\t"+writeHits+"\t"+(normalized ? "exact_modal_support" : "none")+"\t"+sketchMode+
			(legacyK5 ? "" : "\t"+geometry.name())+
			(registerBits==32 ? "\t2\t32\t"+exponent+"\t"+DDLFamilyProfiles.ENCODING32+"\t"+
				(indexed ? "none" : counted ? "scalar_weighted" : Shared.SIMD ? "simd32" : "scalar32") : "")+
			DDLFamilyProfiles.positionMetadataValues(positionMode,partitions,overlap,geometry)+
			(blacklist==null ? "" : "\t"+blacklist.sha80()+"\t"+blacklist.size()+"\t"+blacklist.metadataValues())+"\n";
		Files.write(output.resolve("summary.tsv"),summary.getBytes(StandardCharsets.US_ASCII));
		final StringBuilder histogram=new StringBuilder((normalized ? "floor_support_fraction_times_width" : counted ? "floor_frequency_score_units" : "matching_registers")+"\ttop_hit_queries\texpected_family_queries\n");
		for(int score=0; score<=width; score++){
			histogram.append(score).append('\t').append(topHistogram[score]).append('\t').append(expectedHistogram[score]).append('\n');
		}
		Files.write(output.resolve("score_histogram.tsv"),histogram.toString().getBytes(StandardCharsets.US_ASCII));
		Files.write(output.resolve("recipe.txt"),(String.join("\n",args)+"\n").getBytes(StandardCharsets.UTF_8));
		Files.write(output.resolve("PASS"),"DDL_FAMILY_SEARCH_PASS\n".getBytes(StandardCharsets.US_ASCII));
	}

	private static String required(HashMap<String,String> options,String key){
		final String value=options.get(key);
		if(value==null || value.length()==0){throw new IllegalArgumentException("Missing "+key);}
		return value;
	}

	/** Load benchmark for the four experimental exports, not the production DDLLoader.
	 * Timed sections allocate and retain the same representation; no index or scoring.
	 * Native BGZF is explicitly single-threaded. Equality checks and pins are untimed.
	 */
	private static void loadBenchmark(String[] args) throws Exception{
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg : args){
			final String[] parts=Parse.splitOnFirst(arg,'=');
			if(parts.length!=2 || options.put(parts[0].toLowerCase(Locale.ROOT),parts[1])!=null){throw new IllegalArgumentException("Expected unique flag=value");}
		}
		for(String key : options.keySet()){
			if(!Arrays.asList("in","out","families","seconds","reverse").contains(key)){throw new IllegalArgumentException("Unknown loadbench flag: "+key);}
		}
		final Path input=Paths.get(required(options,"in")),output=Paths.get(required(options,"out"));
		final int families=Integer.parseInt(required(options,"families"));
		final double seconds=Double.parseDouble(options.containsKey("seconds") ? options.get("seconds") : "12");
		if(families<1 || !Double.isFinite(seconds) || seconds<10){throw new IllegalArgumentException("Positive family count and >=10s timed loads required");}
		final boolean reverse=options.containsKey("reverse") && Parse.parseBoolean(options.get("reverse"));
		Shared.setThreads(1);
		Files.createDirectory(output);
		final String[] encodings=reverse ? new String[]{"a48","a48-stream16","a48-triplet","u16"} : new String[]{"u16","a48-triplet","a48-stream16","a48"};
		final LoadedProfiles reference=loadEncoded(input.resolve("library_u16.bgz"),families);
		final StringBuilder table=new StringBuilder("encoding\tfamilies\twidth\tcounts\tloads\tload_seconds\tseconds_per_load\tcompressed_bytes\tinput_sha80\n");
		for(String encoding : encodings){
			final Path path=input.resolve("library_"+encoding+".bgz");
			final String before=sha80(path);
			final LoadedProfiles warm=loadEncoded(path,families);
			if(!warm.encoding.equals(encoding)){throw new IllegalArgumentException("Encoding/header mismatch: "+path);}
			reference.requireEqual(warm);
			long elapsed=0; int loads=0;
			while(elapsed<seconds*1e9){
				final long start=System.nanoTime();
				final LoadedProfiles loaded=loadEncoded(path,families);
				elapsed+=System.nanoTime()-start;
				reference.requireEqual(loaded); loads++;
			}
			pin(path,before);
			final double duration=elapsed*1e-9;
			assert(loads>0 && duration>=10) : "Performance evidence requires >=10 seconds of actual repeated loading";
			table.append(encoding).append('\t').append(families).append('\t').append(reference.width).append('\t')
				.append(reference.counts!=null).append('\t').append(loads).append('\t').append(duration).append('\t')
				.append(duration/loads).append('\t').append(Files.size(path)).append('\t').append(before).append('\n');
			System.err.println("LOAD_COMPLETE encoding="+encoding+" loads="+loads+" timed_seconds="+duration);
		}
		Files.write(output.resolve("summary.tsv"),table.toString().getBytes(StandardCharsets.US_ASCII));
		Files.write(output.resolve("recipe.txt"),(Arrays.toString(args)+"\nNative BgzfInputStream; one worker; warmed repeated loads; equality outside timer; no index.\n").getBytes(StandardCharsets.US_ASCII));
		Files.write(output.resolve("PASS"),"DDL_ENCODING_LOAD_PARITY_PASS\n".getBytes(StandardCharsets.US_ASCII));
	}

	/** Retain identical arrays and metadata for every encoding; counts stay raw here. */
	public static final class LoadedProfiles{
		LoadedProfiles(String encoding_,int width_,int k_,long seed_,String mode_,int families,boolean counted){
			assert(families>0 && width_>0) : "A retained library needs a nonempty roster and fixed-width rows";
			encoding=encoding_; width=width_; k=k_; seed=seed_; mode=mode_;
			names=new String[families]; members=new int[families]; rows=new char[families][]; counts=counted ? new char[families][] : null;
		}
		public void requireEqual(LoadedProfiles other){
			if(other==null || width!=other.width || k!=other.k || seed!=other.seed || !mode.equals(other.mode) ||
				!Arrays.equals(names,other.names) || !Arrays.equals(members,other.members) || (counts==null)!=(other.counts==null)){
				throw new IllegalArgumentException("Decoded library metadata differs");
			}
			for(int f=0; f<rows.length; f++){
				if(!Arrays.equals(rows[f],other.rows[f]) || (counts!=null && !Arrays.equals(counts[f],other.counts[f]))){
					throw new IllegalArgumentException("Decoded registers/counts differ at family "+f);
				}
			}
		}
		public final String encoding,mode;
		public final int width,k;
		public final long seed;
		public final String[] names;
		public final int[] members;
		public final char[][] rows,counts;
	}

	/** Mixed text/binary framing requires fixed-length binary reads, never line splitting
	 * of u16 payloads (which can contain LF, CR and TAB). Same framing for all formats.
	 * The outer buffer uses bulk native reads, also avoiding BgzfInputStream#002's
	 * known scalar read() behavior at empty BGZF blocks.
	 */
	public static LoadedProfiles loadEncoded(Path path,int families) throws Exception{
		if(families<1){throw new IllegalArgumentException("Require expected family count");}
		try(InputStream raw=new BufferedInputStream(Files.newInputStream(path),65536);
			InputStream in=new BufferedInputStream(new BgzfInputStream(raw),65536)){
			final ByteBuilder line=new ByteBuilder();
			if(!readFramedLine(in,line,4096)){throw new EOFException("Missing library header");}
			final String[] header=new String(line.array,0,line.length,StandardCharsets.US_ASCII).split("\t",-1);
			if(header.length!=9 || !header[0].equals("#DDLFamilyProfiles") || !header[1].equals("version=1") ||
				!header[2].startsWith("encoding=") || !header[3].startsWith("width=") || !header[4].equals("registerbits=16") ||
				!header[5].startsWith("k=") || !header[6].startsWith("seed=") || !header[7].startsWith("mode=")){
				throw new IOException("Invalid experimental library header");
			}
			final String encoding=header[2].substring(9),mode=header[7].substring(5);
			final int width=Integer.parseInt(header[3].substring(6)),k=Integer.parseInt(header[5].substring(2));
			final long seed=Long.parseLong(header[6].substring(5));
			final boolean counted=header[8].equals("counts=raw-a48-scaled-from-members");
			if(!Arrays.asList("u16","a48-triplet","a48-stream16","a48").contains(encoding) || width<1 || width>DDLFamilyProfiles.MAX_BUCKETS ||
				k<1 || k>12 || (!mode.equals("literal") && !mode.equals("nonempty")) || (!counted && !header[8].equals("counts=none"))){
				throw new IOException("Unsupported experimental library recipe");
			}
			final LoadedProfiles result=new LoadedProfiles(encoding,width,k,seed,mode,families,counted);
			final byte[] binary=encoding.equals("u16") ? new byte[width*2] : null;
			for(int f=0; f<families; f++){
				if(!readFramedLine(in,line,4096) || line.length<19){throw new EOFException("Missing family header "+f);}
				final byte[] a=line.array;
				int split=8; while(split<line.length && a[split]!='\t'){split++;}
				if(a[0]!='#' || a[1]!='f' || a[2]!='a' || a[3]!='m' || a[4]!='i' || a[5]!='l' || a[6]!='y' || a[7]!='\t' ||
					split==8 || split+9>=line.length || a[split+1]!='m' || a[split+2]!='e' || a[split+3]!='m' || a[split+4]!='b' ||
					a[split+5]!='e' || a[split+6]!='r' || a[split+7]!='s' || a[split+8]!='='){throw new IOException("Invalid family header "+f);}
				result.names[f]=new String(a,8,split-8,StandardCharsets.UTF_8);
				long members=0;
				for(int i=split+9; i<line.length; i++){
					if(a[i]<'0' || a[i]>'9' || (members=members*10+a[i]-'0')>Integer.MAX_VALUE){throw new IOException("Invalid member divisor");}
				}
				if(members<1){throw new IOException("Nonpositive member divisor");}
				result.members[f]=(int)members;
				final char[] row=new char[width]; result.rows[f]=row;
				if(binary!=null){
					int offset=0;
					while(offset<binary.length){final int n=in.read(binary,offset,binary.length-offset); if(n<1){throw new EOFException("Truncated u16 row");} offset+=n;}
					if(in.read()!='\n'){throw new IOException("Missing u16 row terminator");}
					decodeRegisters(binary,binary.length,encoding,row);
				}else{
					if(!readFramedLine(in,line,width*4)){throw new EOFException("Missing text registers");}
					decodeRegisters(line.array,line.length,encoding,row);
				}
				if(counted){
					if(!readFramedLine(in,line,width*4)){throw new EOFException("Missing raw counts");}
					final char[] counts=new char[width]; result.counts[f]=counts;
					decodeRegisters(line.array,line.length,"a48",counts);
				}
			}
			if(in.read()!=-1){throw new IOException("Trailing records/bytes beyond expected family count");}
			return result;
		}
	}

	private static boolean readFramedLine(InputStream in,ByteBuilder line,int max) throws IOException{
		assert(max>0) : "Bounded text framing must not absorb a corrupt binary payload";
		line.clear();
		for(int value=in.read(); value!=-1; value=in.read()){
			if(value=='\n'){return true;}
			if(line.length>=max || value=='\r'){throw new IOException("Invalid or oversized framed line");}
			line.append((byte)value);
		}
		if(line.length!=0){throw new EOFException("Unterminated framed line");}
		return false;
	}

	/** Strict decoder; variable A48 has TAB separators; fixed formats have none. */
	public static void decodeRegisters(byte[] input,int length,String encoding,char[] row) throws IOException{
		if(input==null || row==null || row.length<1 || length<0 || length>input.length){throw new IOException("Invalid decoder buffer");}
		int cursor=0;
		if(encoding.equals("u16")){
			if(length!=row.length*2){throw new IOException("Wrong u16 byte count");}
			for(int b=0; b<row.length; b++){row[b]=(char)((input[cursor++]&255)|((input[cursor++]&255)<<8));}
		}else if(encoding.equals("a48-triplet")){
			if(length!=row.length*3){throw new IOException("Wrong triplet byte count");}
			for(int b=0; b<row.length; b++){
				final int high=symbol(input[cursor++]);
				if(high>15){throw new IOException("Triplet exceeds uint16");}
				row[b]=(char)((high<<12)|(symbol(input[cursor++])<<6)|symbol(input[cursor++]));
			}
		}else if(encoding.equals("a48-stream16")){
			if(length!=(row.length*16+5)/6){throw new IOException("Wrong packed A48 byte count");}
			int pending=0,used=0;
			for(int b=0; b<row.length; b++){
				while(used<16){pending=(pending<<6)|symbol(input[cursor++]); used+=6;}
				used-=16; row[b]=(char)(pending>>>used); pending&=(1<<used)-1;
			}
			if(pending!=0){throw new IOException("Nonzero final A48 padding");}
		}else if(encoding.equals("a48")){
			for(int b=0; b<row.length; b++){
				int value=0,digits=0;
				while(cursor<length && input[cursor]!='\t'){
					if(++digits>3){throw new IOException("Oversized A48 integer");}
					value=(value<<6)|symbol(input[cursor++]);
				}
				if(digits==0 || value>65535){throw new IOException("Missing or oversized A48 value");}
				row[b]=(char)value;
				if(b+1<row.length){if(cursor>=length || input[cursor++]!='\t'){throw new IOException("Missing A48 delimiter");}}
			}
		}else{throw new IOException("Unknown encoding: "+encoding);}
		if(cursor!=length){throw new IOException("Trailing register bytes");}
	}
	private static int symbol(byte value) throws IOException{
		final int symbol=value-48;
		if(symbol<0 || symbol>63){throw new IOException("Invalid A48 symbol");}
		return symbol;
	}

	public static Index loadForSearch(Path list,int width,String mode,int k,long seed,boolean indexed,boolean counted,boolean normalized) throws Exception{
		return load(list,width,mode,k,seed,DDLFamilyProfiles.SKETCH_PLAIN,indexed,counted,normalized,null,null);
	}

	public static Index loadForSearch(Path list,int width,String mode,ProteinKmerGeometry geometry,long seed,
			boolean indexed,boolean counted,boolean normalized) throws Exception{
		return load(list,width,mode,geometry,seed,DDLFamilyProfiles.SKETCH_PLAIN,indexed,counted,normalized,null,null,false);
	}

	public static Index loadForSearch(Path list,int width,String mode,int k,long seed,String sketchMode,
		boolean indexed,boolean counted,boolean normalized) throws Exception{
		final ProteinPositionMode positionMode=ProteinPositionMode.fromSketch(sketchMode);
		if(positionMode.positional() &&
			((width!=512 && width!=1024) || k!=DDLFamilyProfiles.POSITIONAL_K || seed!=0)){
			throw new IllegalArgumentException("Positional sketch v1 search requires buckets=512/1024 k=5 seed=0");
		}
		return load(list,width,mode,k,seed,sketchMode,indexed,counted,normalized,null,null);
	}

	/** Diagnostic-only row retention; identical scoring, roster order and tie rules. */
	public static Index loadForAttribution(Path list,int width,int k,long seed) throws Exception{
		return load(list,width,"nonempty",ProteinKmerGeometry.exact(k),seed,DDLFamilyProfiles.SKETCH_PLAIN,true,false,false,null,null,true);
	}

	/** Diagnostic-only row retention; identical scoring, roster order and tie rules. */
	public static Index loadForAttribution(Path list,int width,int k,long seed,String sketchMode) throws Exception{
		return loadForAttribution(list,width,k,seed,sketchMode,4,8);
	}

	/** Retains rows only after checking the same positional recipe used by the query kernel. */
	public static Index loadForAttribution(Path list,int width,int k,long seed,String sketchMode,int partitions,int overlap) throws Exception{
		return loadForAttribution(list,width,k,seed,sketchMode,partitions,overlap,null);
	}

	/** Bind the exact raw-AA exclusion set while retaining rows for collision diagnostics. */
	public static Index loadForAttribution(Path list,int width,int k,long seed,String sketchMode,int partitions,int overlap,DDLPentamerBlacklist blacklist) throws Exception{
		return loadForAttribution(list,width,ProteinKmerGeometry.exact(k),seed,sketchMode,partitions,overlap,blacklist);
	}

	/** Bind a projected protein geometry while retaining rows for collision diagnostics. */
	public static Index loadForAttribution(Path list,int width,ProteinKmerGeometry geometry,long seed,String sketchMode,
			int partitions,int overlap,DDLPentamerBlacklist blacklist) throws Exception{
		ProteinPositionMode.fromSketch(sketchMode);
		return load(list,width,"nonempty",geometry,seed,sketchMode,false,false,false,null,null,true,partitions,overlap,blacklist);
	}

	/** Explicit structural path: track file stability for inherited or unmeasured payload identities. */
	public static Index loadForAttributionStructural(Path list,int width,ProteinKmerGeometry geometry,long seed,String sketchMode,
			int partitions,int overlap,StructuralInputs inputs) throws Exception{
		return loadForAttributionStructural(list,width,geometry,seed,sketchMode,partitions,overlap,null,inputs);
	}

	/** Structural loading retains the same immutable filter as producer and query workers. */
	public static Index loadForAttributionStructural(Path list,int width,ProteinKmerGeometry geometry,long seed,String sketchMode,
			int partitions,int overlap,DDLPentamerBlacklist blacklist,StructuralInputs inputs)throws Exception{
		if(inputs==null){throw new IllegalArgumentException("Structural loading requires an input tracker");}
		ProteinPositionMode.fromSketch(sketchMode);
		return load(list,width,"nonempty",geometry,seed,sketchMode,false,false,false,null,null,true,partitions,overlap,blacklist,inputs);
	}

	/** Wide schema is mandatory even for an all-zero library; no width inference from values. */
	public static Index loadForSearch32(Path list,int width,String mode,int k,long seed,int exponent,
		boolean indexed,boolean counted,boolean normalized) throws Exception{
		DDLFamilyProfiles.validateLayout(32,exponent);
		if(list==null || width<1 || width>DDLFamilyProfiles.MAX_BUCKETS || Integer.bitCount(width)!=1 || k<1 || k>12 ||
			(!mode.equals("literal") && !mode.equals("nonempty")) || (normalized && (!counted || !mode.equals("nonempty")))){
			throw new IllegalArgumentException("Invalid wide profile search recipe");
		}
		final ArrayList<int[]> rows=new ArrayList<int[]>(),raw=normalized ? new ArrayList<int[]>() : null;
		final ArrayList<char[]> weights=counted && !normalized ? new ArrayList<char[]>() : null;
		final ArrayList<String> names=new ArrayList<String>();
		final HashSet<String> paths=new HashSet<String>(),unique=new HashSet<String>();
		final LineParser1 parser=new LineParser1('\t');
		final int[] support=new int[width],occupancy=new int[width];
		final String kText=Integer.toString(k),seedText=Long.toString(seed),exponentText=Integer.toString(exponent);
		long bytes=0;
		for(String entry : Files.readAllLines(list,StandardCharsets.UTF_8)){
			final String[] fields=entry.split("\t",-1);
			if(fields.length!=3){throw new IllegalArgumentException("Require profile_path, profile_sha80, summary_sha80");}
			Path path=Paths.get(fields[0]);
			if(!path.isAbsolute()){path=list.toAbsolutePath().getParent().resolve(path);}
			path=path.normalize();
			if(!paths.add(path.toString())){throw new IllegalArgumentException("Duplicate wide profile path");}
			pin(path,fields[1]);
			final Path summary=path.getParent().resolve("summary.tsv"); pin(summary,fields[2]);
			final java.util.List<String> metadata=Files.readAllLines(summary,StandardCharsets.US_ASCII);
			if(metadata.size()!=2){throw new IllegalArgumentException("Wide summary requires header and one data row");}
			final String[] keys=metadata.get(0).split("\t",-1),values=metadata.get(1).split("\t",-1);
			if(keys.length!=values.length){throw new IllegalArgumentException("Wide summary field count");}
			final HashMap<String,String> meta=new HashMap<String,String>();
			for(int i=0; i<keys.length; i++){
				if(meta.put(keys[i],values[i])!=null){throw new IllegalArgumentException("Duplicate wide summary field: "+keys[i]);}
			}
			final String namespace=meta.get("namespace");
			validateBlacklistMetadata(meta,null);
			if(namespace==null || !namespace.matches("[A-Za-z0-9_-]+") || !"2".equals(meta.get("version")) ||
				!"32".equals(meta.get("registerbits")) || !DDLFamilyProfiles.ENCODING32.equals(meta.get("encoding")) ||
				!kText.equals(meta.get("k")) || !seedText.equals(meta.get("seed")) || !exponentText.equals(meta.get("exponent"))){
				throw new IllegalArgumentException("Wide summary schema/recipe mismatch");
			}
			final ByteFile file=ByteFile.makeByteFile(path.toString(),false);
			try{
				final byte[] header=file.nextLine();
				if(header==null || !DDLFamilyProfiles.HEADER32.equals(new String(header,StandardCharsets.US_ASCII))){throw new IllegalArgumentException("Wide profile header/schema mismatch");}
				for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
					parser.set(line);
					if(parser.terms()!=16 || !parser.termEquals(namespace,0) || !parser.termEquals("2",10) || !parser.termEquals("32",11) ||
						!parser.termEquals(kText,12) || !parser.termEquals(seedText,13) || !parser.termEquals(exponentText,14) ||
						!parser.termEquals(DDLFamilyProfiles.ENCODING32,15)){
						throw new IllegalArgumentException("Mixed or missing wide row schema/recipe");
					}
					final int rowWidth=parser.parseInt(2),members=parser.parseInt(3);
					if(rowWidth<1 || rowWidth>DDLFamilyProfiles.MAX_BUCKETS || Integer.bitCount(rowWidth)!=1 || members<1 ||
						(!parser.termEquals("literal",6) && !parser.termEquals("nonempty",6))){throw new IllegalArgumentException("Invalid wide row dimensions/mode");}
					if(rowWidth!=width || !parser.termEquals(mode,6)){continue;}
					final String id=parser.parseString(1),name=namespace+":"+id;
					if(id.isEmpty() || !unique.add(name)){throw new IllegalArgumentException("Empty/duplicate wide profile: "+name);}
					final int[] row=new int[width];
					parseUnsignedArray(parser,7,row); parseUnsignedArray(parser,8,support); parseUnsignedArray(parser,9,occupancy);
					final char[] weight=weights==null ? null : new char[width];
					for(int b=0; b<width; b++){
						final int value=row[b],votes=support[b],occupied=occupancy[b];
						if((value>>>(32-exponent))>Math.min(64,(1<<exponent)-1) || votes<0 || votes>members || occupied<0 || occupied>members ||
							(value!=0 && (votes<1 || votes>occupied)) || (value==0 && votes!=(mode.equals("literal") ? members-occupied : 0)) ||
							(mode.equals("literal") && value!=0 && votes<=members-occupied) ||
							(mode.equals("nonempty") && (value==0)!=(occupied==0))){throw new IllegalArgumentException("Invalid wide register/support/occupancy at bucket "+b);}
						if(weight!=null){weight[b]=value==0 ? 0 : decodeFrequency(encodeCount(votes,members),members);}
					}
					rows.add(row); names.add(name);
					if(raw!=null){raw.add(support.clone());}
					if(weights!=null){weights.add(weight);}
				}
			}finally{if(file.close()){throw new IllegalStateException("Wide profile read failed: "+path);}}
			pin(path,fields[1]); pin(summary,fields[2]); bytes+=Files.size(path);
		}
		if(rows.isEmpty()){throw new IllegalArgumentException("No matching wide profiles");}
		final Index result=new Index(names.toArray(new String[names.size()]),rows.toArray(new int[rows.size()][]),
			weights==null ? null : weights.toArray(new char[weights.size()][]),raw==null ? null : raw.toArray(new int[raw.size()][]),indexed);
		// No binary export is emitted: zero packedBytes means unmeasured, not a hypothetical format size.
		result.gzipBytes=bytes;
		return result;
	}

	/** Decimal arrays are parsed directly from bytes, with strict uint32 range before narrowing. */
	private static void parseUnsignedArray(LineParser1 parser,int term,int[] out){
		assert(parser!=null && out!=null) : "A retained primitive row is required before parsing its comma-separated field";
		parser.setBounds(term); final byte[] line=parser.line(); final int end=parser.b();
		int slot=0,digits=0; long value=0;
		for(int p=parser.a(); p<=end; p++){
			final int c=p==end ? ',' : line[p];
			if(c==','){
				if(digits==0 || slot==out.length){throw new IllegalArgumentException("Wide array field count/empty value");}
				out[slot++]=(int)value; value=0; digits=0;
			}else{
				if(c<'0' || c>'9' || digits>=10){throw new IllegalArgumentException("Unsigned decimal32 required");}
				value=value*10+c-'0'; digits++;
				if(value>0xffffffffL){throw new IllegalArgumentException("Unsigned32 overflow");}
			}
		}
		if(slot!=out.length){throw new IllegalArgumentException("Wide array width mismatch");}
	}

	private static Index load(Path list,int width,String mode,int k,long seed,boolean indexed,boolean counted,boolean normalized,Path export,Path formats) throws Exception{
		return load(list,width,mode,ProteinKmerGeometry.exact(k),seed,DDLFamilyProfiles.SKETCH_PLAIN,indexed,counted,normalized,export,formats,false);
	}

	private static Index load(Path list,int width,String mode,int k,long seed,String sketchMode,
		boolean indexed,boolean counted,boolean normalized,Path export,Path formats) throws Exception{
		return load(list,width,mode,ProteinKmerGeometry.exact(k),seed,sketchMode,indexed,counted,normalized,export,formats,false);
	}

	private static Index load(Path list,int width,String mode,ProteinKmerGeometry geometry,long seed,String sketchMode,
		boolean indexed,boolean counted,boolean normalized,Path export,Path formats,boolean retainRows) throws Exception{
		return load(list,width,mode,geometry,seed,sketchMode,indexed,counted,normalized,export,formats,retainRows,4,8);
	}

	private static Index load(Path list,int width,String mode,ProteinKmerGeometry geometry,long seed,String sketchMode,
		boolean indexed,boolean counted,boolean normalized,Path export,Path formats,boolean retainRows,int partitions,int overlap) throws Exception{
		return load(list,width,mode,geometry,seed,sketchMode,indexed,counted,normalized,export,formats,retainRows,partitions,overlap,null);
	}

	/** Shared loader checks exclusion identity before accepting any profile row. */
	private static Index load(Path list,int width,String mode,ProteinKmerGeometry geometry,long seed,String sketchMode,
		boolean indexed,boolean counted,boolean normalized,Path export,Path formats,boolean retainRows,int partitions,int overlap,DDLPentamerBlacklist blacklist) throws Exception{
		return load(list,width,mode,geometry,seed,sketchMode,indexed,counted,normalized,export,formats,retainRows,partitions,overlap,blacklist,null);
	}

	/** Shared loader checks exclusion identity before accepting any profile row. */
	private static Index load(Path list,int width,String mode,ProteinKmerGeometry geometry,long seed,String sketchMode,
		boolean indexed,boolean counted,boolean normalized,Path export,Path formats,boolean retainRows,int partitions,int overlap,
		DDLPentamerBlacklist blacklist,StructuralInputs structural) throws Exception{
		final int k=geometry.k();
		blacklist=ProteinBlacklistPolicy.resolve(geometry,blacklist);
		if(structural!=null && (export!=null || formats!=null)){
			throw new IllegalArgumentException("Structural attribution loading requires metadata-preserving profiles and does not support profile exports");
		}
		if(blacklist!=null && (geometry.span()<5 || export!=null || formats!=null)){
			throw new IllegalArgumentException("Raw pentamer exclusion requires a five-residue span and metadata-preserving profile storage");
		}
		if(!geometry.isExact() && (export!=null || formats!=null)){
			throw new IllegalArgumentException("Irregular amino-acid profiles require TSV metadata");
		}
		final boolean positional=ProteinPositionMode.fromSketch(sketchMode).positional();
		if(positional){
			DDLFamilyProfiles.validatePositionalLayout(partitions,overlap,geometry);
			if((width!=512 && width!=1024) || seed!=0 || export!=null || formats!=null){
				throw new IllegalArgumentException("Positional loading requires width512/1024, seed0 and metadata-preserving input");
			}
		}
		final ArrayList<char[]> rows=new ArrayList<char[]>();
		final ArrayList<char[]> weights=counted && !normalized ? new ArrayList<char[]>() : null;
		final ArrayList<char[]> positions=new ArrayList<char[]>();
		final ArrayList<int[]> rawSupport=normalized ? new ArrayList<int[]>() : null;
		final ArrayList<char[]> diskCounts=counted ? new ArrayList<char[]>() : null;
		final IntList memberCounts=new IntList();
		final ArrayList<String> names=new ArrayList<String>();
		final HashSet<String> unique=new HashSet<String>();
		final HashSet<String> paths=new HashSet<String>();
		if(structural!=null){structural.track(list);}
		for(String entry : Files.readAllLines(list,StandardCharsets.UTF_8)){
			final String[] fields=entry.split("\t",-1);
			if(fields.length!=3 || !paths.add(fields[0])){throw new IllegalArgumentException("Require unique profile_path<TAB>profile_sha80<TAB>summary_sha80 rows");}
			Path path=Paths.get(fields[0]);
			if(!path.isAbsolute()){path=list.getParent().resolve(path);}
			if(structural==null){pin(path,fields[1]);}
			else{requireStructuralIdentity(fields[1]); structural.track(path);}
			final Path summary=path.getParent().resolve("summary.tsv");
			if(structural==null){pin(summary,fields[2]);}
			else{requireStructuralIdentity(fields[2]); structural.track(summary);}
			final java.util.List<String> metadata=Files.readAllLines(summary,StandardCharsets.US_ASCII);
			if(metadata.size()!=2){throw new IllegalArgumentException("Profile summary requires header and one data row");}
			final String[] namesIn=metadata.get(0).split("\t",-1), valuesIn=metadata.get(1).split("\t",-1);
			if(namesIn.length!=valuesIn.length){throw new IllegalArgumentException("Malformed summary");}
			final HashMap<String,String> meta=new HashMap<String,String>();
			for(int i=0; i<namesIn.length; i++){
				if(meta.put(namesIn[i],valuesIn[i])!=null){throw new IllegalArgumentException("Duplicate summary field: "+namesIn[i]);}
			}
			if(meta.containsKey("registerbits") && !"16".equals(meta.get("registerbits"))){throw new IllegalArgumentException("Legacy reader requires registerbits=16");}
			if(!Integer.toString(k).equals(meta.get("k")) || !Long.toString(seed).equals(meta.get("seed")) ||
				!"6".equals(meta.get("exponent"))){throw new IllegalArgumentException("Profile/query sketch recipe mismatch");}
			validateGeometryMetadata(meta,geometry);
			final String profileSketch=meta.containsKey("sketch") ? meta.get("sketch") : DDLFamilyProfiles.SKETCH_PLAIN;
			if(!profileSketch.equals(sketchMode)){throw new IllegalArgumentException("Profile/query sketch mode mismatch");}
			validatePositionalMetadata(meta,sketchMode,partitions,overlap,geometry);
			validateBlacklistMetadata(meta,blacklist);
			final ByteFile file=ByteFile.makeByteFile(path.toString(),false);
			Boolean positioned=null;
			try{
				final byte[] header=file.nextLine();
				if(header==null || !new String(header,StandardCharsets.US_ASCII).startsWith("namespace\tfamily_id\tbuckets\t")){
					throw new IllegalArgumentException("Missing profile header: "+path);
				}
				for(byte[] bytes=file.nextLine(); bytes!=null; bytes=file.nextLine()){
					final String[] parts=new String(bytes,StandardCharsets.US_ASCII).split("\t",-1);
					if(parts.length!=10 && parts.length!=11){throw new IllegalArgumentException("Profile field count");}
					if(positioned==null){positioned=parts.length==11;}
					else if(positioned.booleanValue()!=(parts.length==11)){throw new IllegalArgumentException("Mixed profile position fields");}
					if(!parts[0].equals(meta.get("namespace"))){throw new IllegalArgumentException("Profile namespace mismatch");}
					final int rowWidth=Integer.parseInt(parts[2]);
					if(positional && rowWidth!=width){
						throw new IllegalArgumentException("Positional profile width must match the query width: "+rowWidth+" != "+width);
					}
					// Historical plain-profile tables can contain several independently built widths.
					if(rowWidth!=width || !parts[6].equals(mode)){continue;}
					final String name=parts[0]+":"+parts[1];
					if(!unique.add(name)){throw new IllegalArgumentException("Duplicate profile: "+name);}
					final String[] values=parts[7].split(",",-1);
					if(values.length!=width){throw new IllegalArgumentException("Wrong profile width");}
					final char[] row=new char[width];
					final char[] position=positioned.booleanValue() ? new char[width] : null;
					final char[] weight=counted && !normalized ? new char[width] : null, stored=counted ? new char[width] : null;
					final int[] raw=normalized ? new int[width] : null;
					final String[] support=counted ? parts[8].split(",",-1) : null;
					final String[] positionFields=position==null ? null : parts[10].split(",",-1);
					final int members=Integer.parseInt(parts[3]);
					if(members<1 || (counted && support.length!=width) || (positionFields!=null && positionFields.length!=width)){
						throw new IllegalArgumentException("Invalid member/support/position dimensions");
					}
					for(int b=0; b<width; b++){
						final int value=Integer.parseInt(values[b]);
						if(value<0 || value>65535){throw new IllegalArgumentException("Unsigned16 register required");}
						row[b]=(char)value;
						if(position!=null){
							final int p=Integer.parseInt(positionFields[b]);
							if(p<0 || p>65535 || (value==0)!=(p==0)){throw new IllegalArgumentException("Invalid position at bucket "+b);}
							position[b]=(char)p;
						}
						if(counted){
							final int votes=Integer.parseInt(support[b]);
							if(votes<0 || votes>members){throw new IllegalArgumentException("Support must be within member count");}
							stored[b]=encodeCount(votes,members);
							if(normalized){raw[b]=votes;}
							else{weight[b]=value==0 ? 0 : decodeFrequency(stored[b],members);}
						}
					}
					rows.add(row); names.add(name);
					memberCounts.add(members);
					if(position!=null){positions.add(position);}
					if(counted){if(normalized){rawSupport.add(raw);}else{weights.add(weight);} diskCounts.add(stored);}
				}
			}finally{if(file.close()){throw new IllegalStateException("Failed reading profile file: "+path);}}
			if(structural==null){
				pin(path,fields[1]); pin(summary,fields[2]);
			}else{
				structural.track(path); structural.track(summary);
			}
		}
		if(rows.isEmpty()){throw new IllegalArgumentException("No matching profiles");}
		final String[] nameArray=names.toArray(new String[names.size()]);
		final char[][] rowArray=rows.toArray(new char[rows.size()][]);
		final long[] sizes=exportProfiles(export,nameArray,rowArray,diskCounts,memberCounts,k,seed,mode);
		if(formats!=null){exportFormats(formats,nameArray,rowArray,diskCounts,memberCounts,k,seed,mode);}
		final Index result=normalized ? new Index(nameArray,rowArray,rawSupport.toArray(new int[rawSupport.size()][]),indexed,true) :
			new Index(nameArray,rowArray,counted ? weights.toArray(new char[weights.size()][]) : null,
				positions.isEmpty() ? null : positions.toArray(new char[positions.size()][]),indexed);
		result.packedBytes=sizes[0]; result.gzipBytes=sizes[1];
		if(retainRows){result.attributionRows=rowArray;}
		if(structural!=null){structural.verify();}
		return result;
	}

	/** Missing fields denote the historical unfiltered recipe; partial or mismatched fields fail. */
	public static void validateBlacklistMetadata(java.util.Map<String,String> meta,DDLPentamerBlacklist blacklist){
		if(meta==null){throw new IllegalArgumentException("Profile metadata is required");}
		final String pin=meta.get("blacklist_sha80"), words=meta.get("blacklist_words"),
			k=meta.get("blacklist_k"), mode=meta.get("blacklist_mode"), provenance=meta.get("blacklist_provenance");
		if(blacklist==null){
			if(pin!=null || words!=null || k!=null || mode!=null || provenance!=null){
				throw new IllegalArgumentException("Blacklisted profiles require the same query blacklist");
			}
		}else if(!blacklist.sha80().equals(pin) || !Integer.toString(blacklist.size()).equals(words)){
			throw new IllegalArgumentException("Profile/query pentamer blacklist identity or size mismatch");
		}else{
			int present=0;
			if(k!=null){present++;}
			if(mode!=null){present++;}
			if(provenance!=null){present++;}
			if(present==0 && blacklist.wordLength()==5 && blacklist.matchMode()==DDLPentamerBlacklist.MatchMode.PREFIX &&
					!blacklist.sha80().equals("unmeasured")){
				return;
			}
			if(present!=3 || !Integer.toString(blacklist.wordLength()).equals(k) ||
					!blacklist.matchMode().label().equals(mode) || !blacklist.provenance().equals(provenance)){
				throw new IllegalArgumentException("Profile/query protein blacklist mode or provenance mismatch");
			}
		}
	}

	/** Legacy positional summaries encode only the original4/8 recipe; partial metadata is invalid. */
	public static void validatePositionalMetadata(java.util.Map<String,String> meta,String sketchMode,int partitions,int overlap){
		validatePositionalMetadata(meta,sketchMode,partitions,overlap,ProteinKmerGeometry.exact(DDLFamilyProfiles.POSITIONAL_K));
	}

	public static void validatePositionalMetadata(java.util.Map<String,String> meta,String sketchMode,int partitions,int overlap,ProteinKmerGeometry geometry){
		if(meta==null){throw new IllegalArgumentException("Profile metadata is required to bind query sketch layout");}
		final ProteinPositionMode mode=ProteinPositionMode.fromSketch(sketchMode);
		int present=0;
		for(String key:new String[]{"partitions","overlap","tag_bits","segment_words"}){
			if(meta.containsKey(key)){present++;}
		}
		final String basis=meta.get("position_basis");
		if(!mode.positional()){
			if(present!=0 || basis!=null){throw new IllegalArgumentException("Plain profiles cannot declare positional layout fields");}
			return;
		}
		DDLFamilyProfiles.validatePositionalLayout(partitions,overlap,geometry,mode);
		if(mode==ProteinPositionMode.RELATIVE){
			if(present!=3 || meta.containsKey("segment_words") || !DDLFamilyProfiles.RELATIVE_POSITION_BASIS.equals(basis) ||
					!Integer.toString(partitions).equals(meta.get("partitions")) || !Integer.toString(overlap).equals(meta.get("overlap")) ||
					!Integer.toString(DDLFamilyProfiles.positionalTagBits(partitions)).equals(meta.get("tag_bits"))){
				throw new IllegalArgumentException("Relative profile/query layout requires complete residue-length metadata: "+meta);
			}
			return;
		}
		if(basis!=null){throw new IllegalArgumentException("Fixed profiles cannot declare a relative position basis");}
		if(present==0 && mode==ProteinPositionMode.FIXED){
			if(!geometry.isExact(DDLFamilyProfiles.POSITIONAL_K) || partitions!=4 || overlap!=8){
				throw new IllegalArgumentException("Legacy positional metadata supports only exact K5 partitions=4 overlap=8");
			}
			return;
		}
		if(present!=4 || !Integer.toString(partitions).equals(meta.get("partitions")) ||
				!Integer.toString(overlap).equals(meta.get("overlap")) ||
				!Integer.toString(DDLFamilyProfiles.positionalTagBits(partitions)).equals(meta.get("tag_bits")) ||
				!Integer.toString(DDLFamilyProfiles.positionalSegmentWords(partitions,geometry)).equals(meta.get("segment_words"))){
			throw new IllegalArgumentException("Profile/query positional layout mismatch or incomplete metadata: "+meta);
		}
	}

	public static void validateGeometryMetadata(java.util.Map<String,String> meta,ProteinKmerGeometry geometry){
		if(meta==null || geometry==null){throw new IllegalArgumentException("Profile geometry metadata is required");}
		final String stored=meta.get("geometry");
		if(stored==null){
			if(!geometry.isExact(5)){throw new IllegalArgumentException("Missing profile amino-acid geometry metadata");}
			return;
		}
		if(!geometry.name().equals(stored)){throw new IllegalArgumentException("Profile/query amino-acid geometry mismatch");}
	}

	/** Compare only register encodings: identical metadata and raw-count layout in all files. */
	private static void exportFormats(Path out,String[] names,char[][] rows,ArrayList<char[]> counts,
		IntList members,int k,long seed,String mode) throws Exception{
		assert(rows.length>0 && names.length==rows.length) : "Every serialized profile must have a named family";
		final boolean force=ReadWrite.FORCE_BGZIP,nativeOut=ReadWrite.PREFER_NATIVE_BGZF_OUT;
		final boolean allow=ReadWrite.ALLOW_NATIVE_BGZF,bgzf=ReadWrite.USE_BGZF,external=ReadWrite.USE_BGZIP;
		final int zip=ReadWrite.ZIPLEVEL;
		final StringBuilder stats=new StringBuilder("encoding\twidth\tfamilies\tcounts\tcompression\traw_bytes\tcompressed_bytes\n");
		try{
			ReadWrite.FORCE_BGZIP=ReadWrite.PREFER_NATIVE_BGZF_OUT=ReadWrite.ALLOW_NATIVE_BGZF=ReadWrite.USE_BGZF=true;
			ReadWrite.USE_BGZIP=false; ReadWrite.ZIPLEVEL=5;
			for(String encoding : new String[]{"u16","a48-triplet","a48-stream16","a48"}){
				final Path path=out.resolve("library_"+encoding+".bgz");
				final ByteStreamWriter writer=new ByteStreamWriter(path.toString(),false,false,false);
				final ByteBuilder buffer=new ByteBuilder();
				buffer.append("#DDLFamilyProfiles\tversion=1\tencoding=").append(encoding).append("\twidth=").append(rows[0].length)
					.append("\tregisterbits=16\tk=").append(k).append("\tseed=").append(seed).append("\tmode=").append(mode)
					.append("\tcounts=").append(counts==null ? "none" : "raw-a48-scaled-from-members").nl();
				long bytes=buffer.length;
				writer.start();
				try{
					writer.print(buffer);
					for(int f=0; f<rows.length; f++){
						buffer.clear(); buffer.append("#family\t").append(names[f]).append("\tmembers=").append(members.get(f)).nl();
						encodeRegisters(rows[f],encoding,buffer); buffer.nl();
						if(counts!=null){encodeRegisters(counts.get(f),"a48",buffer); buffer.nl();}
						bytes+=buffer.length; writer.print(buffer);
					}
				}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Packed format export failed");}}
				stats.append(encoding).append('\t').append(rows[0].length).append('\t').append(rows.length).append('\t')
					.append(counts!=null).append("\tnative-bgzf-level5\t").append(bytes).append('\t').append(Files.size(path)).append('\n');
			}
		}finally{
			ReadWrite.FORCE_BGZIP=force; ReadWrite.PREFER_NATIVE_BGZF_OUT=nativeOut;
			ReadWrite.ALLOW_NATIVE_BGZF=allow; ReadWrite.USE_BGZF=bgzf; ReadWrite.USE_BGZIP=external; ReadWrite.ZIPLEVEL=zip;
		}
		Files.write(out.resolve("formats.tsv"),stats.toString().getBytes(StandardCharsets.US_ASCII));
	}

	/** A48 triplets pad each16-bit key to18bits; stream16 packs BEFORE A48 encoding.
	 * The stream is MSB-first, padded only at its final6-bit character. Binary baselines
	 * remain little endian and are explicitly distinct from the requested A48 formats.
	 */
	public static void encodeRegisters(char[] row,String encoding,ByteBuilder buffer){
		assert(row!=null && buffer!=null) : "Serialization requires complete register values and an output buffer";
		if(encoding.equals("a48")){
			for(int b=0; b<row.length; b++){if(b>0){buffer.tab();} buffer.appendA48((int)row[b]);}
		}else if(encoding.equals("a48-triplet")){
			for(char value : row){
				buffer.append((byte)(48+(value>>>12))).append((byte)(48+((value>>>6)&63))).append((byte)(48+(value&63)));
			}
		}else if(encoding.equals("a48-stream16")){
			long pending=0; int used=0;
			for(char value : row){
				pending=(pending<<16)|value; used+=16;
				while(used>=6){used-=6; buffer.append((byte)(48+((pending>>>used)&63)));}
				pending&=(1L<<used)-1;
			}
			if(used>0){buffer.append((byte)(48+((pending<<(6-used))&63)));}
		}else{
			final int bits=encoding.equals("u16") ? 16 : encoding.equals("packed24") ? 24 : encoding.equals("packed18") ? 18 : -1;
			if(bits<0){throw new IllegalArgumentException("Unknown register encoding: "+encoding);}
			long pending=0; int used=0;
			for(char value : row){
				pending|=(long)value<<used; used+=bits;
				while(used>=8){buffer.append((byte)pending); pending>>>=8; used-=8;}
			}
			if(used>0){buffer.append((byte)pending);}
		}
	}

	/** Raw uint16 on disk; scale is derivable from the preserved true divisor. */
	public static char encodeCount(int votes,int members){
		if(members<1 || votes<0 || votes>members){throw new IllegalArgumentException("Invalid support/member count");}
		final long divisor=((long)members+65534)/65535;
		final long code=((long)votes+divisor/2)/divisor;
		assert(code>=0 && code<=65535) : "Scaling by ceil(members/65535) bounds every raw support code";
		return (char)Math.min(65535,code);
	}

	/** Load-time Q16 frequency, including reconstruction of scaled large-family counts. */
	public static char decodeFrequency(char code,int members){
		if(members<1){throw new IllegalArgumentException("Positive family divisor required");}
		final long divisor=((long)members+65534)/65535;
		final long frequency=((long)code*divisor*65535+members/2)/members;
		assert(frequency>=0) : "Count reconstruction and Q16 multiplication must remain in long arithmetic";
		return (char)Math.min(65535,frequency);
	}

	/** Experimental compact download payload; little endian, fixed register widths.
	 * Header: magic8, flags/width/k/seedLow/seedHigh/families uint32. Each row:
	 * UTF8-name byte length uint32, name bytes, optional true members uint32,
	 * unsigned16 DDL values, optional unsigned16 RAW scaled counts (not frequencies).
	 */
	private static long[] exportProfiles(Path output,String[] names,char[][] rows,ArrayList<char[]> counts,
		IntList members,int k,long seed,String mode) throws Exception{
		assert(names.length==rows.length && names.length==members.size && rows.length>0) : "Serialized names, divisors and profiles must have the same roster";
		final ByteStreamWriter writer=output==null ? null : new ByteStreamWriter(output.toString(),false,false,false);
		final ByteBuilder buffer=new ByteBuilder();
		buffer.append("DDFP0001").appendUint32((counts==null ? 0 : 1)|(mode.equals("nonempty") ? 2 : 0))
			.appendUint32(rows[0].length).appendUint32(k).appendUint32(seed&0xffffffffL).appendUint32(seed>>>32).appendUint32(rows.length);
		long bytes=buffer.length;
		if(writer!=null){writer.start(); writer.print(buffer);}
		for(int f=0; f<rows.length; f++){
			final byte[] name=names[f].getBytes(StandardCharsets.UTF_8);
			bytes+=4+name.length+(counts==null ? 2L : 4L)*rows[f].length+(counts==null ? 0 : 4);
			if(writer==null){continue;}
			buffer.clear(); buffer.appendUint32(name.length).append(name);
			if(counts!=null){buffer.appendUint32(members.get(f));}
			for(char value : rows[f]){buffer.append((byte)value).append((byte)(value>>>8));}
			if(counts!=null){for(char value : counts.get(f)){buffer.append((byte)value).append((byte)(value>>>8));}}
			writer.print(buffer);
		}
		if(writer!=null && writer.poisonAndWait()){throw new IllegalStateException("Failed writing compact profile payload");}
		return new long[]{bytes,output==null ? -1 : Files.size(output)};
	}

	private static void pin(Path path,String expected) throws Exception{
		requirePin(expected);
		if(!sha80(path).equals(expected)){throw new IllegalStateException("Input sha80 mismatch: "+path);}
	}
	private static void requirePin(String expected){
		if(expected==null || !expected.matches("[0-9a-f]{20}")){throw new IllegalArgumentException("Require sha80 pin");}
	}
	/** Only explicit structural loading may accept an identity that was not measured. */
	private static void requireStructuralIdentity(String expected){
		if(!"unmeasured".equals(expected)){requirePin(expected);}
	}
	private static String sha80(Path path) throws Exception{
		final MessageDigest md=MessageDigest.getInstance("SHA-256");
		final byte[] buffer=new byte[65536];
		try(InputStream stream=Files.newInputStream(path)){
			for(int n=stream.read(buffer); n!=-1; n=stream.read(buffer)){md.update(buffer,0,n);}
		}
		final byte[] digest=md.digest();
		final StringBuilder actual=new StringBuilder(20);
		for(int i=digest.length-10; i<digest.length; i++){
			final int value=digest[i]&255;
			actual.append("0123456789abcdef".charAt(value>>>4)).append("0123456789abcdef".charAt(value&15));
		}
		return actual.toString();
	}

	/** Tracks observable file stability; inherited pins or unmeasured markers do not certify bytes. */
	public static final class StructuralInputs{
		public Path track(Path path) throws IOException{
			final Path key=normalize(path);
			final FileState old=tracked.get(key);
			if(old!=null){
				verify(key,old);
			}else{
				tracked.put(key,FileState.snapshot(key));
			}
			return key;
		}

		public void verify() throws IOException{
			for(java.util.Map.Entry<Path,FileState> entry:tracked.entrySet()){
				verify(entry.getKey(),entry.getValue());
			}
		}

		private static void verify(Path path,FileState expected) throws IOException{
			final FileState actual=FileState.snapshot(path);
			if(actual.size!=expected.size || !actual.modified.equals(expected.modified) || !Objects.equals(actual.key,expected.key)){
				throw new IllegalStateException("Structural input changed while loading: "+path);
			}
		}

		private static Path normalize(Path path){
			if(path==null){throw new IllegalArgumentException("Structural input path is required");}
			return path.toAbsolutePath().normalize();
		}

		private final HashMap<Path,FileState> tracked=new HashMap<Path,FileState>();
	}

	private static final class FileState{
		static FileState snapshot(Path path) throws IOException{
			final BasicFileAttributes attrs=Files.readAttributes(path,BasicFileAttributes.class,LinkOption.NOFOLLOW_LINKS);
			if(!attrs.isRegularFile()){throw new IllegalArgumentException("Require regular non-symlink input file: "+path);}
			return new FileState(attrs.size(),attrs.lastModifiedTime(),attrs.fileKey());
		}
		FileState(long size_,Object modified_,Object key_){
			size=size_; modified=modified_; key=key_;
		}
		final long size;
		final Object modified,key;
	}

	/** Dense key offsets plus contiguous family postings: no object per key/hit. */
	public static final class Index{
		/** Builds the uint16 CSR by disjoint bucket ranges, preserving ascending family postings.
		 * Inputs are borrowed only until this method returns; callers must not mutate them during construction.
		 * Existing constructors remain serial reference implementations and support the other index modes.
		 */
		public static Index parallel(String[] names,char[][] profiles,char[][] weights,int threads){
			return parallel(names,profiles,weights,null,threads);
		}

		public static Index parallel(String[] names,char[][] profiles,char[][] weights,char[][] positions,int threads){
			final long start=System.nanoTime();
			if(names==null || profiles==null || names.length==0 || names.length!=profiles.length ||
					profiles[0]==null || threads<1 || (weights!=null && weights.length!=profiles.length) ||
					(positions!=null && positions.length!=profiles.length)){
				throw new IllegalArgumentException("Parallel CSR requires named uint16 rows and positive thread count");
			}
			final int width=profiles[0].length;
			if(width<1 || width>DDLFamilyProfiles.MAX_BUCKETS || Integer.bitCount(width)!=1){
				throw new IllegalArgumentException("Invalid parallel CSR width: "+width);
			}
			for(int f=0; f<profiles.length; f++){
				if(profiles[f]==null || profiles[f].length!=width ||
						(weights!=null && (weights[f]==null || weights[f].length!=width)) ||
						(positions!=null && (positions[f]==null || positions[f].length!=width))){
					throw new IllegalArgumentException("Mixed parallel CSR widths at family "+f);
				}
			}
			final int workers=Math.min(threads,width);
			final int[] offsets=new int[(width<<16)+1];
			final BuildShard[] shards=new BuildShard[workers];
			for(int i=0; i<workers; i++){
				shards[i]=new BuildShard(profiles,weights,positions,offsets,
					(int)((long)width*i/workers),(int)((long)width*(i+1)/workers));
			}
			final long allocated=System.nanoTime();
			final ExecutorService pool=workers==1 ? null : Executors.newFixedThreadPool(workers);
			try{
				runBuildShards(pool,shards,false,null,null,null,null);
				final long counted=System.nanoTime();
				long count=0, nonemptyKeys=0; int empty=0;
				for(BuildShard shard:shards){
					shard.base=(int)count; count+=shard.count;
					nonemptyKeys+=shard.nonemptyKeys;
					if(count>Integer.MAX_VALUE-8){throw new IllegalArgumentException("Posting count exceeds array capacity");}
				}
				for(int f=0; f<profiles.length; f++){
					boolean present=false;
					for(BuildShard shard:shards){present|=shard.nonempty[f];}
					if(!present){empty++;}
				}
				final int[] postings=new int[(int)count],cursor=new int[offsets.length-1];
				final char[] postingWeights=weights==null ? null : new char[(int)count];
				final char[] postingPositions=positions==null ? null : new char[(int)count];
				final long prepared=System.nanoTime();
				runBuildShards(pool,shards,true,cursor,postings,postingWeights,postingPositions);
				final long filled=System.nanoTime();
				assert(offsets[offsets.length-1]==postings.length) : "Disjoint CSR ranges must cover every nonzero positive-weight register";
				long countWork=0,prefixWork=0,cursorWork=0,fillWork=0;
				for(BuildShard shard:shards){
					countWork+=shard.countNanos; prefixWork+=shard.prefixNanos;
					cursorWork+=shard.cursorNanos; fillWork+=shard.fillNanos;
				}
				final BuildTiming timing=new BuildTiming(workers,allocated-start,counted-allocated,
					prepared-counted,filled-prepared,countWork,prefixWork,cursorWork,fillWork);
				return new Index(names,width,empty,offsets,postings,postingWeights,postingPositions,nonemptyKeys,timing);
			}finally{if(pool!=null){pool.shutdown();}}
		}

		/** Optional direct-completion rows; the ordinary builder retains none.
		 * The short rows are owned copies, so callers may release their build inputs.
		 */
		public static Index parallel(String[] names,char[][] profiles,char[][] weights,int threads,boolean retainDirectRows){
			if(retainDirectRows && weights!=null){throw new IllegalArgumentException("Direct completion requires unweighted profiles");}
			final Index index=parallel(names,profiles,weights,threads);
			if(retainDirectRows){
				index.directRows=new short[profiles.length][index.width];
				for(int family=0; family<profiles.length; family++){
					for(int bin=0; bin<index.width; bin++){index.directRows[family][bin]=(short)profiles[family][bin];}
				}
			}
			return index;
		}

		/** Inspect optional allocation without materializing a second row set. */
		public boolean hasDirectRows(){return directRows!=null;}

		/** Compare caller-owned rows to this index's real unweighted uint16 CSR. */
		public boolean matchesProfiles(char[][] rows){
			if(rows==null || rows.length!=names.length || postingWeights!=null || postingPositions!=null ||
					offsets==null || postings==null || offsets.length!=(width<<16)+1){
				return false;
			}
			final int[] expected=new int[Character.MAX_VALUE+1];
			long total=0;
			for(int b=0; b<width; b++){
				Arrays.fill(expected, 0);
				for(char[] row:rows){
					if(row==null || row.length!=width){return false;}
					final int value=row[b];
					if(value>0){expected[value]++; total++;}
				}
				final int base=b<<16;
				if(offsets[base]!=offsets[base+1]){return false;}
				for(int value=1; value<expected.length; value++){
					final int start=offsets[base|value], end=offsets[(base|value)+1];
					if(end<start || end-start!=expected[value]){return false;}
					int last=-1;
					for(int p=start; p<end; p++){
						final int family=postings[p];
						if(family<=last || family<0 || family>=rows.length || rows[family][b]!=(char)value){return false;}
						last=family;
					}
				}
			}
			return total==postings.length;
		}

		/** Construct one explicitly experimental policy before publishing the index to workers. */
		public static Index parallelPruned(String[] names,char[][] profiles,int threads,boolean early,boolean direct,String order){
			if(early && direct){throw new IllegalArgumentException("Early rejection and direct completion are separate experimental arms");}
			final int orderCode=parseScanOrder(order);
			final Index index=parallel(names,profiles,null,threads,direct);
			index.pruned=early || direct; index.directCompletion=direct; index.scanOrder=orderCode;
			return index;
		}

		/** Stable option encoding:0=interleaved,1=start quadrants,2=stop quadrants,3=plain bins. */
		public static int parseScanOrder(String order){
			if("plain".equals(order)){return 3;}
			if("interleaved".equals(order)){return 0;}
			if("startquadrant".equals(order)){return 1;}
			if("stopquadrant".equals(order)){return 2;}
			throw new IllegalArgumentException("indexorder must be plain, interleaved, startquadrant or stopquadrant");
		}
		public boolean usesPruning(){return pruned;}
		public boolean usesDirectCompletion(){return directCompletion;}
		/** Whether the ordinary score(char[], Workspace) path can accumulate inverse
		 * posting depth alongside raw hits in the same complete CSR traversal.
		 * An empty index has this capability but has mean depth zero; bindings that
		 * require a measured positive normalization must reject that separately. */
		public boolean supportsInversePostingRarity(){
			return !pruned && !directCompletion && wideRows==null && offsets!=null && postings!=null &&
				postingWeights==null && totalSupport==null;
		}
		public int scanOrder(){return scanOrder;}

		/** Every submitted task is joined even after a peer fails or the caller is interrupted. */
		private static void runBuildShards(ExecutorService pool,BuildShard[] shards,boolean fill,
				int[] cursor,int[] postings,char[] weights,char[] positions){
			assert(shards.length>0) : "A nonempty bucket partition must have at least one worker";
			if(pool==null){
				assert(shards.length==1) : "Only the one-worker path executes on the caller thread";
				if(fill){shards[0].fill(cursor,postings,weights,positions);}else{shards[0].count();}
				return;
			}
			final ArrayList<Future<?>> jobs=new ArrayList<>(shards.length);
			Throwable failure=null; boolean interrupted=false;
			try{
				for(final BuildShard shard:shards){
					jobs.add(pool.submit(()->{if(fill){shard.fill(cursor,postings,weights,positions);}else{shard.count();}}));
				}
			}catch(Throwable problem){failure=problem;}
			for(Future<?> job:jobs){
				boolean done=false;
				while(!done){
					try{job.get(); done=true;}
					catch(InterruptedException problem){interrupted=true; if(failure==null){failure=problem;}}
					catch(ExecutionException problem){done=true; if(failure==null){failure=problem.getCause();}}
				}
			}
			if(interrupted){Thread.currentThread().interrupt();}
			if(failure instanceof Error){throw (Error)failure;}
			if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
			if(failure!=null){throw new IllegalStateException("Parallel CSR construction failed",failure);}
		}

		private Index(String[] names,int width,int empty,int[] offsets,int[] postings,char[] weights,char[] positions,
				long nonemptyKeys,BuildTiming timing){
			assert(offsets.length==(width<<16)+1 && offsets[offsets.length-1]==postings.length) : "Published CSR must retain its complete key and posting ranges";
			this.names=names.clone(); this.width=width; emptyProfiles=empty;
			this.offsets=offsets; this.postings=postings; postingWeights=weights;
			postingPositions=positions; references=null; referenceWeights=null; referencePositions=null; buildTiming=timing;
			meanPostingDepth=DDLFamilySearch.meanPostingDepth(postings.length, nonemptyKeys);
		}

		/** Wall phases are sequential; worker sums overlap and must not be added to wall time. */
		public static final class BuildTiming{
			private BuildTiming(int threads,long setup,long countPrefix,long allocate,long fill,
					long countWork,long prefixWork,long cursorWork,long fillWork){
				this.threads=threads; setupWallNanos=setup; countPrefixWallNanos=countPrefix;
				allocationWallNanos=allocate; fillWallNanos=fill; countWorkerNanos=countWork;
				prefixWorkerNanos=prefixWork; cursorWorkerNanos=cursorWork; fillWorkerNanos=fillWork;
			}
			public final int threads;
			public final long setupWallNanos,countPrefixWallNanos,allocationWallNanos,fillWallNanos;
			public final long countWorkerNanos,prefixWorkerNanos,cursorWorkerNanos,fillWorkerNanos;
		}

		/** Only this shard writes counts at (firstKey,lastKey] and cursors at [firstKey,lastKey). */
		private static final class BuildShard{
			BuildShard(char[][] rows,char[][] weights,char[][] positions,int[] offsets,int firstBucket,int lastBucket){
				assert(firstBucket>=0 && lastBucket>firstBucket && lastBucket<=rows[0].length) : "Bucket partitions must be nonempty and disjoint";
				this.rows=rows; this.weights=weights; this.positions=positions; this.offsets=offsets;
				first=firstBucket; last=lastBucket; nonempty=new boolean[rows.length];
			}
			void count(){
				final long begin=System.nanoTime();
				for(int f=0; f<rows.length; f++){
					final char[] row=rows[f],weight=weights==null ? null : weights[f];
					boolean present=false;
					for(int b=first; b<last; b++){
						if(row[b]>0 && (weight==null || weight[b]>0)){
							offsets[((b<<16)|row[b])+1]++; count++; present=true;
						}
					}
					nonempty[f]=present;
				}
				final long counted=System.nanoTime();
				if(count>Integer.MAX_VALUE-8){throw new IllegalArgumentException("CSR shard exceeds posting capacity");}
				int sum=0;
				for(int i=(first<<16)+1,end=last<<16; i<=end; i++){
					final int n=offsets[i]; nonemptyKeys+=(n>0 ? 1 : 0); sum+=n; offsets[i]=sum;
				}
				assert(sum==count) : "Local prefix sum must conserve this shard's register count";
				countNanos=counted-begin; prefixNanos=System.nanoTime()-counted;
			}
			void fill(int[] cursor,int[] postings,char[] postingWeights,char[] postingPositions){
				final long begin=System.nanoTime();
				final int firstKey=first<<16,lastKey=last<<16;
				cursor[firstKey]=base;
				for(int i=firstKey+1; i<lastKey; i++){cursor[i]=base+offsets[i];}
				for(int i=firstKey+1; i<=lastKey; i++){offsets[i]+=base;}
				final long prepared=System.nanoTime();
				final int firstPosting=base,lastPosting=(int)(base+count);
				for(int f=0; f<rows.length; f++){
					final char[] row=rows[f],weight=weights==null ? null : weights[f];
					final char[] position=positions==null ? null : positions[f];
					for(int b=first; b<last; b++){
						if(row[b]>0 && (weight==null || weight[b]>0)){
							final int p=cursor[(b<<16)|row[b]]++;
							assert(p>=firstPosting && p<lastPosting) : "A bucket owner must fill only its own posting slice";
							postings[p]=f;
							if(weight!=null){postingWeights[p]=weight[b];}
							if(position!=null){postingPositions[p]=position[b];}
						}
					}
				}
				cursorNanos=prepared-begin; fillNanos=System.nanoTime()-prepared;
			}
			final char[][] rows,weights,positions;
			final int[] offsets;
			final int first,last;
			final boolean[] nonempty;
			long nonemptyKeys;
			int base;
			long count,countNanos,prefixNanos,cursorNanos,fillNanos;
		}

		/** Wide rows are immutable copies; prefix postings are candidates, never matches. */
		public Index(String[] names,int[][] profiles,char[][] weights,int[][] support,boolean indexed){
			buildTiming=null;
			if(names==null || profiles==null || names.length==0 || names.length!=profiles.length || profiles[0]==null){
				throw new IllegalArgumentException("Nonempty named wide profiles required");
			}
			width=profiles[0].length;
			if(width<1 || width>DDLFamilyProfiles.MAX_BUCKETS || Integer.bitCount(width)!=1 ||
				(weights!=null && weights.length!=names.length) || (support!=null && support.length!=names.length) || (weights!=null && support!=null)){
				throw new IllegalArgumentException("Wide dimensions or mutually exclusive weight/support mismatch");
			}
			this.names=names.clone(); wideRows=new int[profiles.length][];
			wideWeights=weights==null ? null : new char[weights.length][];
			wideSupport=support==null ? null : new int[support.length][];
			totalSupport=support==null ? null : new long[support.length];
			offsets=indexed ? new int[(width<<16)+1] : null;
			references=null; referenceWeights=null; referencePositions=null; postingWeights=null; postingPositions=null;
			long count=0; int empty=0;
			for(int f=0; f<profiles.length; f++){
				if(profiles[f]==null || profiles[f].length!=width || (weights!=null && (weights[f]==null || weights[f].length!=width)) ||
					(support!=null && (support[f]==null || support[f].length!=width))){throw new IllegalArgumentException("Mixed wide profile dimensions");}
				wideRows[f]=profiles[f].clone();
				if(weights!=null){wideWeights[f]=weights[f].clone();}
				if(support!=null){wideSupport[f]=support[f].clone();}
				final long before=count;
				for(int b=0; b<width; b++){
					final int value=profiles[f][b];
					if(support!=null){
						final int n=support[f][b];
						if(n<0 || (value==0)!=(n==0)){throw new IllegalArgumentException("Nonempty wide modal support must be positive, empty support zero");}
						totalSupport[f]+=n;
					}
					if(value!=0 && (weights==null || weights[f][b]!=0)){
						if(indexed){offsets[((b<<16)|(value>>>16))+1]++;} count++;
					}
				}
				if(count==before){empty++;}
			}
			emptyProfiles=empty;
			if(!indexed){postings=null; meanPostingDepth=0; return;}
			if(count>Integer.MAX_VALUE-8){throw new IllegalArgumentException("Wide posting count exceeds array capacity");}
			long nonemptyKeys=0;
			for(int i=1; i<offsets.length; i++){
				final int n=offsets[i]; nonemptyKeys+=(n>0 ? 1 : 0); offsets[i]=n+offsets[i-1];
			}
			postings=new int[(int)count]; final int[] cursor=offsets.clone();
			for(int f=0; f<wideRows.length; f++){
				for(int b=0; b<width; b++){
					final int value=wideRows[f][b];
					if(value!=0 && (wideWeights==null || wideWeights[f][b]!=0)){postings[cursor[(b<<16)|(value>>>16)]++]=f;}
				}
			}
			assert(offsets[offsets.length-1]==postings.length) : "Wide CSR covers all positive-weight nonempty registers, including prefix0";
			meanPostingDepth=DDLFamilySearch.meanPostingDepth(count, nonemptyKeys);
		}

		public void score(int[] query,Workspace work){
			if(wideRows==null || query==null || query.length!=width || work.scores.length!=names.length){throw new IllegalArgumentException("Wide query/workspace dimensions or registerbits mismatch");}
			if(work.rarityScores!=null || work.lastRarityScores!=null){throw new IllegalArgumentException("Posting-rarity scoring is only supported by the raw uint16 CSR path");}
			if(totalSupport!=null && work.rawScores==null){work.rawScores=new long[names.length];}
			for(int i=0; i<work.size; i++){
				final int f=work.touched[i]; work.scores[f]=0; if(work.rawScores!=null){work.rawScores[f]=0;}
			}
			work.size=0; work.postings=0; work.comparedRegisters=0;
			if(offsets==null){
				for(int f=0; f<names.length; f++){
					long sum=0;
					if(Shared.SIMD && wideSupport==null && wideWeights==null){sum=SIMDFamilyProfiles.matches(query,wideRows[f]);}
					else{
						for(int b=0; b<width; b++){
							if(query[b]!=0 && query[b]==wideRows[f][b]){sum+=wideSupport!=null ? wideSupport[f][b] : wideWeights!=null ? wideWeights[f][b] : 1;}
						}
					}
					if(sum>0){
						if(totalSupport!=null){work.rawScores[f]=sum;}else{work.scores[f]=(int)sum;}
						work.touched[work.size++]=f;
					}
				}
				work.comparedRegisters=(long)width*names.length;
			}else{
				for(int b=0; b<width; b++){
					final int value=query[b]; if(value==0){continue;}
					final int key=(b<<16)|(value>>>16),end=offsets[key+1];
					work.postings+=end-offsets[key];
					for(int p=offsets[key]; p<end; p++){
						final int f=postings[p];
						if(value!=wideRows[f][b]){continue;}
						if(totalSupport!=null){
							if(work.rawScores[f]==0){work.touched[work.size++]=f;}
							work.rawScores[f]+=wideSupport[f][b];
						}else{
							if(work.scores[f]==0){work.touched[work.size++]=f;}
							work.scores[f]+=wideWeights==null ? 1 : wideWeights[f][b];
						}
					}
				}
			}
			assert(work.size<=names.length) : "Only exact nonzero uint32 matches produce unique touched families";
		}

		public int[] profileCopy32(int family){
			if(wideRows==null){throw new IllegalStateException("A legacy profile cannot be widened into a member-built mode");}
			return wideRows[family].clone();
		}

		public Index(String[] names,char[][] profiles){
			this(names,profiles,null,true);
		}

		/** Exact raw modal-support fractions, with no Q16 or uint16 count rounding.
		 * All nonzero modes must have positive support; empty modes have zero support.
		 * The scan oracle uses scalar long sums; existing SIMD paths are unchanged.
		 */
		public Index(String[] names,char[][] profiles,int[][] support,boolean indexed,boolean normalized){
			this(names,profiles,null,indexed);
			if(!normalized || support==null || support.length!=profiles.length){throw new IllegalArgumentException("Require one raw support row per family");}
			totalSupport=new long[profiles.length];
			postingSupport=indexed ? new int[postings.length] : null;
			referenceSupport=indexed ? null : new int[profiles.length][];
			final int[] cursor=indexed ? offsets.clone() : null;
			for(int f=0; f<profiles.length; f++){
				if(support[f]==null || support[f].length!=width){throw new IllegalArgumentException("Raw support width mismatch");}
				for(int b=0; b<width; b++){
					final int count=support[f][b],value=profiles[f][b];
					if(count<0 || (value==0)!=(count==0)){throw new IllegalArgumentException("Nonempty modal support must be positive, empty support zero");}
					totalSupport[f]+=count;
					if(indexed && value!=0){
						final int p=cursor[(b<<16)|value]++;
						assert(postings[p]==f) : "Raw supports must use the identical family/bucket CSR ordering";
						postingSupport[p]=count;
					}
				}
				if(!indexed){referenceSupport[f]=support[f].clone();}
			}
		}

		/** The two backends use identical positive scores and ordinal tie breaks. */
		public Index(String[] names,char[][] profiles,char[][] weights,boolean indexed){
			this(names,profiles,weights,null,indexed);
		}

		public Index(String[] names,char[][] profiles,char[][] weights,char[][] positions,boolean indexed){
			buildTiming=null;
			if(names.length==0 || names.length!=profiles.length){throw new IllegalArgumentException("Nonempty named profiles required");}
			width=profiles[0].length;
			if(width<1 || width>DDLFamilyProfiles.MAX_BUCKETS || Integer.bitCount(width)!=1){throw new IllegalArgumentException("Invalid width");}
			if((weights!=null && weights.length!=profiles.length) || (positions!=null && positions.length!=profiles.length)){
				throw new IllegalArgumentException("Weight/position roster mismatch");
			}
			this.names=names.clone();
			offsets=indexed ? new int[(width<<16)+1] : null;
			references=indexed ? null : new short[profiles.length][];
			referenceWeights=indexed || weights==null ? null : new short[profiles.length][];
			referencePositions=indexed || positions==null ? null : new short[profiles.length][];
			long count=0;
			int empty=0;
			for(int f=0; f<profiles.length; f++){
				final char[] row=profiles[f], weight=weights==null ? null : weights[f];
				final char[] position=positions==null ? null : positions[f];
				if(row.length!=width || (weight!=null && weight.length!=width) || (position!=null && position.length!=width)){
					throw new IllegalArgumentException("Mixed widths");
				}
				final long before=count;
				for(int b=0; b<width; b++){
					if(row[b]>0 && (weight==null || weight[b]>0)){
						if(indexed){offsets[((b<<16)|row[b])+1]++;}
						count++;
					}
				}
				if(!indexed){
					references[f]=shorts(row);
					if(weight!=null){referenceWeights[f]=shorts(weight);}
					if(position!=null){referencePositions[f]=shorts(position);}
				}
				if(count==before){empty++;}
			}
			emptyProfiles=empty;
			if(!indexed){postings=null; postingWeights=null; postingPositions=null; meanPostingDepth=0; return;}
			if(count>Integer.MAX_VALUE-8){throw new IllegalArgumentException("Posting count exceeds array capacity");}
			long nonemptyKeys=0;
			for(int i=1; i<offsets.length; i++){
				final int n=offsets[i]; nonemptyKeys+=(n>0 ? 1 : 0); offsets[i]=n+offsets[i-1];
			}
			postings=new int[(int)count];
			postingWeights=weights==null ? null : new char[(int)count];
			postingPositions=positions==null ? null : new char[(int)count];
			final int[] cursor=offsets.clone();
			for(int f=0; f<profiles.length; f++){
				final char[] row=profiles[f], weight=weights==null ? null : weights[f];
				final char[] position=positions==null ? null : positions[f];
				for(int b=0; b<width; b++){
					if(row[b]>0 && (weight==null || weight[b]>0)){
						final int p=cursor[(b<<16)|row[b]]++;
						postings[p]=f;
						if(weight!=null){postingWeights[p]=weight[b];}
						if(position!=null){postingPositions[p]=position[b];}
					}
				}
			}
			assert(offsets[offsets.length-1]==postings.length) : "CSR offsets must cover each positive-weight nonzero register once";
			meanPostingDepth=DDLFamilySearch.meanPostingDepth(count, nonemptyKeys);
		}

		public void score(char[] query,Workspace work){
			if(wideRows!=null || query.length!=width || work.scores.length!=names.length){throw new IllegalArgumentException("Query/workspace dimensions or registerbits mismatch");}
			final double[] rarityScores=work.rarityScores;
			if(rarityScores!=null && (rarityScores.length!=names.length || offsets==null || postingWeights!=null || totalSupport!=null)){
				throw new IllegalArgumentException("Posting-rarity scoring requires the unweighted indexed uint16 CSR path");
			}
			work.clearStaleRarity();
			if(totalSupport!=null){scoreNormalized(query,work); return;}
			for(int i=0; i<work.size; i++){
				final int family=work.touched[i]; work.scores[family]=0;
				if(work.rawScores!=null){work.rawScores[family]=0;}
				if(rarityScores!=null){rarityScores[family]=0;}
			}
			work.size=0; work.postings=0; work.comparedRegisters=0;
			if(offsets==null){
				if(work.queryShorts==null || work.queryShorts.length!=width){work.queryShorts=new short[width];}
				for(int b=0; b<width; b++){work.queryShorts[b]=(short)query[b];}
				for(int f=0; f<names.length; f++){
					final short[] weight=referenceWeights==null ? null : referenceWeights[f];
					final int score=Shared.SIMD ? (weight==null ? SIMDFamilyProfiles.matches(work.queryShorts,references[f]) :
						SIMDFamilyProfiles.weighted(work.queryShorts,references[f],weight)) : scalar(work.queryShorts,references[f],weight);
					if(score>0){work.scores[f]=score; work.touched[work.size++]=f;}
				}
				work.comparedRegisters=(long)width*names.length;
				assert(work.size<=names.length) : "A full scan must add each positive-scoring family at most once";
				return;
			}
			if(rarityScores==null){
				for(int b=0; b<width; b++){
					if(query[b]==0){continue;}
					final int key=(b<<16)|query[b], end=offsets[key+1];
					work.postings+=end-offsets[key];
					for(int p=offsets[key]; p<end; p++){
						final int family=postings[p];
						if(work.scores[family]==0){work.touched[work.size++]=family;}
						work.scores[family]+=postingWeights==null ? 1 : postingWeights[p];
					}
				}
			}else{
				for(int b=0; b<width; b++){
					if(query[b]==0){continue;}
					final int key=(b<<16)|query[b], start=offsets[key], end=offsets[key+1], count=end-start;
					work.postings+=count;
					if(count==0){continue;}
					final double rarity=meanPostingDepth/count;
					for(int p=start; p<end; p++){
						final int family=postings[p];
						if(work.scores[family]==0){work.touched[work.size++]=family;}
						work.scores[family]++;
						rarityScores[family]+=rarity;
					}
				}
			}
			assert(work.size<=names.length) : "Each candidate enters the touched list once";
		}

		/** Diagnostic-only raw-hit scan in even-occupied-key then odd-occupied-key order.
		 * Completes every posting list and returns the same scores as score(); no early exit.
		 * The ordinary scanner above has no instrumentation branch or progress allocation.
		 */
		public void scoreInterleaved(char[] query,Workspace work,DDLIndexScanProgress progress,int truthFamily){
			scoreCheckpoints(query,work,progress,truthFamily,null);
		}

		/** Diagnostic explicit occupied-bin order; null retains the original even/odd permutation. */
		public void scoreCheckpoints(char[] query,Workspace work,DDLIndexScanProgress progress,int truthFamily,int[] order){
			scoreCheckpoints(query,work,progress,truthFamily,order,null);
		}

		/** Optional explicit boundaries measure whole-quadrant groups, including empty groups. */
		public void scoreCheckpoints(char[] query,Workspace work,DDLIndexScanProgress progress,int truthFamily,int[] order,int[] boundaries){
			if(query==null || work==null || progress==null || query.length!=width || work.scores.length!=names.length ||
					work.rarityScores!=null || work.lastRarityScores!=null || offsets==null || wideRows!=null ||
					totalSupport!=null || postingWeights!=null || truthFamily< -1 || truthFamily>=names.length){
				throw new IllegalArgumentException("Checkpoints require matching uint16 unweighted CSR/query/workspace and optional valid truth family");
			}
			for(int i=0; i<work.size; i++){
				final int family=work.touched[i]; work.scores[family]=0;
				if(work.rawScores!=null){work.rawScores[family]=0;}
			}
			work.size=0; work.postings=0; work.comparedRegisters=0;
			if(order==null){progress.reset(query);}else{progress.reset(query,order);}
			if(boundaries!=null){progress.setBoundaries(boundaries);}
			int best=0,leader=-1,twoPlus=0,consumed=0;
			progress.capture(work,0,best,leader,twoPlus,truthFamily);
			for(int visited=0; visited<progress.keys; visited++){
				final int evens=(progress.keys+1)/2;
				final int i=order!=null ? visited : visited<evens ? 2*visited : 2*(visited-evens)+1;
				final int bin=progress.bins[i],key=(bin<<16)|query[bin],end=offsets[key+1];
				work.postings+=end-offsets[key];
				for(int p=offsets[key]; p<end; p++){
					final int family=postings[p];
					if(work.scores[family]==0){work.touched[work.size++]=family;}
					final int score=++work.scores[family];
					if(score==2){twoPlus++;}
					if(score>best || (score==best && (leader<0 || family<leader))){best=score; leader=family;}
				}
				progress.capture(work,++consumed,best,leader,twoPlus,truthFamily);
			}
			progress.requireComplete();
			assert(consumed==progress.keys && work.size<=names.length && best<=width) :
				"Every occupied query key is visited once; each family contributes at most one hit per register";
		}

		/** Experimental half-scan rejection or direct completion; ordinary score() is unchanged.
		 * order contains each occupied bin once; null scans natural bins and skips zeros.
		 * In the null case keys=width and boundary is a bin boundary, normally width/2.
		 * Early mode returns no candidates when the partial maximum is at most3.
		 * Direct mode instead finishes only families with at least2 partial matches.
		 * @return true only when early mode rejected the entire query
		 */
		public boolean scorePruned(char[] query,Workspace work,int[] order,int keys,int boundary,boolean direct){
			if(query==null || query.length!=width || work==null || work.scores.length!=names.length ||
					keys<0 || keys>width || (order==null ? keys!=width : order.length<keys) || boundary<0 || boundary>keys ||
					work.rarityScores!=null || work.lastRarityScores!=null || offsets==null || postingWeights!=null ||
					totalSupport!=null || wideRows!=null || direct && directRows==null){
				throw new IllegalArgumentException("Pruned search requires raw uint16 CSR, an occupied-bin order and optional retained direct rows");
			}
			for(int i=0; i<work.size; i++){
				final int family=work.touched[i]; work.scores[family]=0;
				if(work.rawScores!=null){work.rawScores[family]=0;}
			}
			work.size=0; work.postings=0; work.comparedRegisters=0;
			int best=0;
			for(int i=0; i<boundary; i++){
				final int bin=order==null ? i : order[i];
				assert(bin>=0 && bin<width && (order==null || query[bin]!=0)) : "A permutation selects occupied bins; a plain scan visits the full bin range";
				if(query[bin]==0){continue;}
				final int key=(bin<<16)|query[bin],end=offsets[key+1];
				work.postings+=end-offsets[key];
				for(int p=offsets[key]; p<end; p++){
					final int family=postings[p];
					if(work.scores[family]==0){work.touched[work.size++]=family;}
					final int score=++work.scores[family];
					if(score>best){best=score;}
				}
			}
			if(direct){
				if(work.queryShorts==null || work.queryShorts.length!=width){work.queryShorts=new short[width];}
				for(int bin=0; bin<width; bin++){work.queryShorts[bin]=(short)query[bin];}
				int survivors=0;
				for(int i=0; i<work.size; i++){
					final int family=work.touched[i];
					if(work.scores[family]<2){work.scores[family]=0; continue;}
					work.scores[family]=Shared.SIMD ? SIMDFamilyProfiles.matches(work.queryShorts,directRows[family]) :
						scalar(work.queryShorts,directRows[family],null);
					assert(work.scores[family]>=2) : "Full matching-register count cannot be lower than the retained partial count";
					work.touched[survivors++]=family;
				}
				work.size=survivors; work.comparedRegisters=(long)survivors*width;
				return false;
			}
			if(best<=3){
				for(int i=0; i<work.size; i++){work.scores[work.touched[i]]=0;}
				work.size=0;
				return true;
			}
			for(int i=boundary; i<keys; i++){
				final int bin=order==null ? i : order[i];
				assert(bin>=0 && bin<width && (order==null || query[bin]!=0)) : "The second half visits the remaining permutation or natural bins";
				if(query[bin]==0){continue;}
				final int key=(bin<<16)|query[bin],end=offsets[key+1];
				work.postings+=end-offsets[key];
				for(int p=offsets[key]; p<end; p++){
					final int family=postings[p];
					if(work.scores[family]==0){work.touched[work.size++]=family;}
					work.scores[family]++;
				}
			}
			assert(work.size<=names.length) : "Each positive-scoring family appears once in the completed shortlist input";
			return false;
		}

		/** Exact single-family scan for family-isolated calibration workers. */
		public int scoreFamily(short[] query,int family){
			if(wideRows!=null || references==null || totalSupport!=null || query==null ||
					query.length!=width || family<0 || family>=names.length){
				throw new IllegalArgumentException("Single-family scoring requires legacy non-indexed raw profiles");
			}
			final short[] weight=referenceWeights==null ? null : referenceWeights[family];
			final int score=Shared.SIMD ? (weight==null ? SIMDFamilyProfiles.matches(query,references[family]) :
				SIMDFamilyProfiles.weighted(query,references[family],weight)) : scalar(query,references[family],weight);
			assert(score>=0 && score<=width*65535) : "Single-family exact score is bounded by the register count and uint16 weights";
			return score;
		}

		/** Descending matches, deterministic family ordinal tie-break, positive scores only. */
		public int top(Workspace work,int[] result){
			if(totalSupport!=null){return topNormalized(work,result);}
			Arrays.fill(result,-1);
			int size=0;
			for(int t=0; t<work.size; t++){
				final int family=work.touched[t], score=work.scores[family];
				// Most touched families cannot enter a full shortlist; avoid scanning
				// all50 ranks just to rediscover that. Preserve the ordinal tie rule.
				if(size>0 && size==result.length){
					final int last=result[size-1];
					if(score<work.scores[last] || (score==work.scores[last] && family>=last)){continue;}
				}
				int rank=0;
				while(rank<size && (work.scores[result[rank]]>score ||
					(work.scores[result[rank]]==score && result[rank]<family))){rank++;}
				if(rank>=result.length){continue;}
				for(int j=Math.min(size,result.length-1); j>rank; j--){result[j]=result[j-1];}
				result[rank]=family;
				if(size<result.length){size++;}
			}
			return size;
		}

		private void scoreNormalized(char[] query,Workspace work){
			assert(totalSupport!=null) : "Normalized scoring requires raw support denominators";
			if(work.rawScores==null){work.rawScores=new long[names.length];}
			for(int i=0; i<work.size; i++){
				final int family=work.touched[i]; work.rawScores[family]=0; work.scores[family]=0;
			}
			work.size=0; work.postings=0; work.comparedRegisters=0;
			if(offsets==null){
				for(int f=0; f<names.length; f++){
					long sum=0;
					for(int b=0; b<width; b++){if(query[b]!=0 && query[b]==(references[f][b]&65535)){sum+=referenceSupport[f][b];}}
					if(sum>0){work.rawScores[f]=sum; work.touched[work.size++]=f;}
				}
				work.comparedRegisters=(long)width*names.length;
			}else{
				for(int b=0; b<width; b++){
					if(query[b]==0){continue;}
					final int key=(b<<16)|query[b],end=offsets[key+1];
					work.postings+=end-offsets[key];
					for(int p=offsets[key]; p<end; p++){
						final int family=postings[p];
						if(work.rawScores[family]==0){work.touched[work.size++]=family;}
						work.rawScores[family]+=postingSupport[p];
					}
				}
			}
			assert(work.size<=names.length) : "Each positive raw numerator must appear once in the candidate list";
		}

		private int topNormalized(Workspace work,int[] result){
			assert(totalSupport!=null && work.rawScores!=null) : "Rank only after raw-support scoring";
			Arrays.fill(result,-1);
			int size=0;
			for(int t=0; t<work.size; t++){
				final int family=work.touched[t];
				if(size==result.length && (size==0 || compareFamilies(work,family,result[size-1])<=0)){continue;}
				int rank=0;
				while(rank<size && compareFamilies(work,result[rank],family)>0){rank++;}
				for(int j=Math.min(size,result.length-1); j>rank; j--){result[j]=result[j-1];}
				result[rank]=family;
				if(size<result.length){size++;}
			}
			return size;
		}

		private int compareFamilies(Workspace work,int a,int b){
			assert(work.rawScores[a]>0 && work.rawScores[b]>0) : "Only positive fractions enter the shortlist";
			final int comparison=compareFractions(work.rawScores[a],totalSupport[a],work.rawScores[b],totalSupport[b]);
			return comparison!=0 ? comparison : Integer.compare(b,a);
		}

		private int histogramBin(Workspace work,int family,int scale){
			if(totalSupport==null){return work.scores[family]/scale;}
			final long numerator=work.rawScores[family],denominator=totalSupport[family];
			assert(numerator>=0 && numerator<=denominator) : "Matched modal support is a subset of total support";
			return denominator==0 ? 0 : (int)(numerator*width/denominator);
		}
		public int families(){return names.length;}
		public double meanPostingDepth(){return meanPostingDepth;}
		/** Read-only census of the built CSR, independent of any query distribution.
		 * Percentiles use nearest rank; nonempty and all-slot denominators are explicit.
		 * No histogram or work is added to construction or ordinary search. */
		public PostingStatistics postingStatistics(){
			if(offsets==null || postings==null){throw new IllegalStateException("Posting statistics require a materialized CSR index");}
			if(offsets.length<2 || offsets[0]!=0 || offsets[offsets.length-1]!=postings.length){
				throw new IllegalStateException("CSR offsets must span exactly the stored posting entries");
			}
			final long[] histogram=new long[names.length+1];
			long entries=0,squares=0,nonempty=0; int maximum=0;
			for(int key=0; key<offsets.length-1; key++){
				final int length=offsets[key+1]-offsets[key];
				if(length<0 || length>names.length){
					throw new IllegalStateException("CSR key "+key+" has invalid length "+length+"; at most one posting per family/key");
				}
				histogram[length]++;
				entries=Math.addExact(entries,length);
				squares=Math.addExact(squares,Math.multiplyExact((long)length,length));
				if(length>0){nonempty++;maximum=Math.max(maximum,length);}
			}
			if(entries!=postings.length){throw new IllegalStateException("CSR census must conserve stored entries");}
			return new PostingStatistics(names.length,offsets.length-1,nonempty,entries,squares,maximum,
				percentile99(histogram,nonempty,1),percentile99(histogram,offsets.length-1,0));
		}
		/** Nearest-rank 99th percentile over the selected histogram population. */
		private static int percentile99(long[] histogram,long count,int first){
			assert(histogram!=null && count>=0 && (first==0 || first==1)) : "CSR percentile populations include either all or only occupied keys";
			if(count==0){return 0;}
			final long rank=(Math.multiplyExact(count,99)+99)/100;
			long cumulative=0;
			for(int length=first; length<histogram.length; length++){
				cumulative=Math.addExact(cumulative,histogram[length]);
				if(cumulative>=rank){return length;}
			}
			throw new IllegalStateException("CSR histogram does not cover percentile population "+count);
		}
		/** Immutable counts; empty populations have mean and percentile zero. */
		public static final class PostingStatistics{
			/** Construct the census after entry conservation and histogram checks. */
			private PostingStatistics(int families_,long keySlots_,long nonemptyLists_,long entries_,long sumSquares_,
					int maximum_,int p99Nonempty_,int p99All_){
				assert(families_>=0 && keySlots_>0 && nonemptyLists_<=keySlots_ && entries_>=nonemptyLists_) :
					"Built CSR statistics must conserve occupied keys and postings";
				families=families_;
				keySlots=keySlots_;
				nonemptyLists=nonemptyLists_;
				entries=entries_;
				sumSquares=sumSquares_;
				maximum=maximum_;
				p99Nonempty=p99Nonempty_;
				p99All=p99All_;
				meanNonempty=nonemptyLists==0 ? 0 : entries/(double)nonemptyLists;
				meanAll=entries/(double)keySlots;
			}
			public final int families,maximum,p99Nonempty,p99All;
			public final long keySlots,nonemptyLists,entries,sumSquares;
			public final double meanNonempty,meanAll;
		}
		/** Caller gets an independent copy, never a mutable scoring reference. */
		public char[] profileCopy(int family){
			if(attributionRows==null){throw new IllegalStateException("Rows were not retained for attribution");}
			return attributionRows[family].clone();
		}
		private static short[] shorts(char[] values){
			assert(values!=null) : "Unsigned profile bit patterns must exist before SIMD conversion";
			final short[] result=new short[values.length];
			for(int i=0; i<values.length; i++){result[i]=(short)values[i];}
			return result;
		}
		public final String[] names;
		public final int emptyProfiles;
		public long packedBytes,gzipBytes;
		private final int width;
		private final double meanPostingDepth;
		private final int[] offsets,postings;
		private final char[] postingWeights,postingPositions;
		private final short[][] references,referenceWeights,referencePositions;
		private char[][] attributionRows;
		public long[] totalSupport;
		private int[] postingSupport;
		private int[][] referenceSupport;
		private int[][] wideRows,wideSupport;
		private char[][] wideWeights;
		/** Null for the historical serial constructors; populated by parallel(). */
		public final BuildTiming buildTiming;
		private short[][] directRows;
		private boolean pruned,directCompletion;
		private int scanOrder;
	}
	private static double meanPostingDepth(long postings,long nonemptyKeys){return nonemptyKeys==0 ? 0 : postings/(double)nonemptyKeys;}

	/** Exact rational ordering without overflowing cross-products or floating ties.
	 * Fast multiplication for ordinary counts; continued fractions for large sums.
	 */
	public static int compareFractions(long a,long b,long c,long d){
		if(a<0 || c<0 || b<1 || d<1){throw new IllegalArgumentException("Nonnegative numerators and positive denominators required");}
		if(a<=Long.MAX_VALUE/d && c<=Long.MAX_VALUE/b){return Long.compare(a*d,c*b);}
		int sign=1;
		while(true){
			final long x=a/b,y=c/d;
			if(x!=y){return sign*Long.compare(x,y);}
			final long r=a%b,s=c%d;
			if(r==0 || s==0){return sign*Long.compare(r,s);}
			a=b; b=r; c=d; d=s; sign=-sign;
		}
	}

	/** Java8-safe fallback; weights are unsigned16, not signed shorts. */
	public static int scalar(short[] query,short[] reference,short[] weights){
		assert(query.length==reference.length && (weights==null || weights.length==query.length) && query.length<=32768) :
			"One reference weight per query bucket; the unsigned sum must fit signed int";
		int score=0;
		for(int b=0; b<query.length; b++){
			if(query[b]!=0 && query[b]==reference[b]){score+=weights==null ? 1 : weights[b]&65535;}
		}
		return score;
	}

	/** Reused query state: clearing is O(touched families), not O(all families). */
	public static final class Workspace{
		public Workspace(int families){scores=new int[families]; touched=new int[families];}
		private void clearStaleRarity(){
			if(lastRarityScores!=null && lastRarityScores!=rarityScores){
				for(int i=0; i<size; i++){lastRarityScores[touched[i]]=0;}
			}
			lastRarityScores=rarityScores;
		}
		public final int[] scores,touched;
		public long[] rawScores;
		public double[] rarityScores;
		private double[] lastRarityScores;
		public int size;
		public long postings,comparedRegisters;
		private short[] queryShorts;
	}
}
