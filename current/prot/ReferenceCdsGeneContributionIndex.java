package prot;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;

import fileIO.ByteFile;
import map.IntHashMap;
import parse.LineParser1;

/**
 * Builds and reads the compact, random-access reference-CDS contribution index.
 *
 * <p>The accepted whole-genome contribution artifact currently has exactly one
 * physical CDS row per gene identity. Schema 1 records that validated condition
 * explicitly and stores each gene ordinal as six bytes: its physical CDS length
 * ({@code int}) followed by its assigned family rank plus one ({@code unsigned
 * short}; zero means no decrementable tracked-family assignment). A future
 * multipart corpus must use a new schema that also retains row count and the sum
 * of squared physical-row lengths; this builder fails rather than flattening it.</p>

 * <p>Contribution rows may interleave TIDs as long as each TID's zero-based gene
 * ordinals remain contiguous. The builder therefore uses two passes: the first
 * validates and sizes each TID's fixed record range, and the second memory-maps
 * the output and writes every {@code (tid, ordinal)} directly into that range.</p>
 *
 * <p>The source manifest and current whole-genome cache SHA-256 values are
 * embedded as raw bytes. The loader also requires caller-supplied hashes for the
 * index, manifest, and cache, preventing a contribution index from being paired
 * with a stale 4,432-family cache or a different gene-call corpus.</p>
 *
 * @author Yoimiya
 */
public final class ReferenceCdsGeneContributionIndex {

	private ReferenceCdsGeneContributionIndex(){}

	/** Builds an index or validates and summarizes an existing one. */
	public static void main(final String[] args) throws Exception{
		final Config c=Config.parse(args);
		if(c.build){
			build(c);
		}else{
			verifySha256(c.manifest,c.manifestSha256,"contribution manifest");
			verifySha256(c.cache,c.cacheSha256,"whole-genome cache");
			try(Reader reader=Reader.load(c.index,c.indexSha256,c.manifestSha256,
					c.cacheSha256,c.expectedFamilies)){
				reportNativeCounts(reader, c);
				System.err.println("ReferenceCdsGeneContributionIndex PASS: tids="+
					reader.tids()+" genes="+reader.genes()+" index="+c.index);
			}
		}
	}

	/**
	 * Reports integer gene-count bounds through the validated native readers.
	 * Optional survival inputs verify identical TID coverage and counts before any
	 * bound is reported. This measures the six-decimal label separation used by
	 * MagQCVectorMaker.appendFmt; it does not classify individual stored vectors.
	 */
	private static void reportNativeCounts(final Reader reader, final Config c){
		assert(reader.counts.length>0) : "Reader.load requires a nonempty gene directory before count summaries";
		final ReferenceCdsShredSurvivalTableReader survival;
		if(c.survivalManifest!=null){
			final String manifestHash=verifySha80(c.survivalManifest, c.survivalSha80);
			final String cacheHash=verifySha80(c.shredsCache, c.shredsSha80);
			survival=ReferenceCdsShredSurvivalTableReader.load(c.survivalManifest, manifestHash, cacheHash);
			if(survival.tidCount()!=reader.tids()){
				throw new IllegalArgumentException("Survival/index TID coverage differs: "+survival.tidCount()+" != "+reader.tids());
			}
		}else{survival=null;}
		int max=0, maxTid=-1;
		long total=0;
		for(int i=0; i<reader.counts.length; i++){
			final int tid=reader.tids[i], count=reader.counts[i];
			if(survival!=null && survival.cdsTotal(tid)!=count){
				throw new IllegalArgumentException("Survival/index native gene count differs for tid="+tid+
					": survival="+survival.cdsTotal(tid)+" index="+count);
			}
			total+=count;
			if(count>max || (count==max && tid<maxTid)){max=count; maxTid=tid;}
		}
		if(total!=reader.genes()){throw new IllegalStateException("Directory count sum differs from validated gene total");}
		// R/N loses at least 1/N when incomplete; F/(R+F)>=1/(N+1) for F>=1, R<=N.
		// Leave the half-unit tie excluded: this is a sufficient, conservative bound.
		final boolean separated=(max<1999999);
		System.err.println("ReferenceCdsGeneContributionIndex counts: tids="+reader.tids()+
			" genes="+total+" max_native_total="+max+" max_tid="+maxTid+
			" min_missing_gene_fraction="+(1.0/max)+" min_positive_contamination="+(1.0/(max+1L))+
			" six_decimal_separated="+separated+" survival_parity="+(survival==null ? "NOT_REQUESTED" : "PASS"));
	}

	/** Checks an operator-facing SHA80 pin and returns the full digest for existing native APIs. */
	private static String verifySha80(final String path, final String expected){
		assert(path!=null && expected!=null) : "Config requires paths and SHA80 pins together for survival count validation";
		final String observed=ReferenceCdsSurvivalLabelReader.sha256File(path);
		if(!observed.endsWith(expected)){throw new IllegalArgumentException("SHA80 mismatch for "+path);}
		return observed;
	}

	/**
	 * Memory-mapped immutable contribution reader. Directory state is small and
	 * heap-resident; the gene records remain in the operating system page cache.
	 */
	public static final class Reader implements AutoCloseable {
		private final FileChannel channel;
		private final MappedByteBuffer data;
		private final long genes;
		private final int[] tids,counts;
		private final long[] firstGenes;
		private final IntHashMap tidToEntry;

		private Reader(final FileChannel channel_, final MappedByteBuffer data_,
				final long genes_, final int[] tids_, final long[] firstGenes_,
				final int[] counts_, final IntHashMap tidToEntry_){
			channel=channel_; data=data_; genes=genes_; tids=tids_;
			firstGenes=firstGenes_; counts=counts_; tidToEntry=tidToEntry_;
		}

		/** Loads and cryptographically binds one schema-1 index. */
		public static Reader load(final String path, final String expectedIndexSha256,
				final String expectedManifestSha256, final String expectedCacheSha256,
				final int expectedFamilies) throws IOException{
			if(expectedFamilies<1 || expectedFamilies>Character.MAX_VALUE){
				throw new IllegalArgumentException("Schema 1 requires expectedFamilies in [1,"+
					(int)Character.MAX_VALUE+"]: "+expectedFamilies);
			}
			requireHex64(expectedIndexSha256,"index SHA-256");
			requireHex64(expectedManifestSha256,"manifest SHA-256");
			requireHex64(expectedCacheSha256,"cache SHA-256");
			verifySha256(path,expectedIndexSha256,"contribution index");
			final FileChannel channel=FileChannel.open(Paths.get(path),StandardOpenOption.READ);
			boolean success=false;
			try{
				final long size=channel.size();
				if(size<HEADER_BYTES || size>Integer.MAX_VALUE){
					throw new IllegalArgumentException("Contribution index size outside schema-1 range: "+size);
				}
				final MappedByteBuffer data=channel.map(FileChannel.MapMode.READ_ONLY,0,size);
				data.order(ByteOrder.BIG_ENDIAN);
				final byte[] magic=new byte[MAGIC.length]; data.get(magic);
				if(!java.util.Arrays.equals(magic,MAGIC)){
					throw new IllegalArgumentException("Contribution index magic mismatch: "+path);
				}
				final int version=data.getInt(),headerBytes=data.getInt();
				final int families=data.getInt(),nTids=data.getInt();
				final long genes=data.getLong(),directoryOffset=data.getLong();
				final byte[] manifestHash=new byte[32],cacheHash=new byte[32];
				data.get(manifestHash); data.get(cacheHash);
				if(version!=VERSION || headerBytes!=HEADER_BYTES || families!=expectedFamilies ||
						nTids<1 || genes<1 || directoryOffset!=HEADER_BYTES+genes*RECORD_BYTES ||
						size!=directoryOffset+(long)nTids*DIRECTORY_BYTES){
					throw new IllegalArgumentException("Contribution index header/count mismatch: "+path);
				}
				if(!expectedManifestSha256.equals(hex(manifestHash)) ||
						!expectedCacheSha256.equals(hex(cacheHash))){
					throw new IllegalArgumentException("Contribution index source binding mismatch: "+path);
				}
				final int[] tids=new int[nTids],counts=new int[nTids];
				final long[] firstGenes=new long[nTids];
				final IntHashMap tidToEntry=new IntHashMap(nTids*2);
				data.position((int)directoryOffset);
				long priorEnd=0;
				for(int i=0; i<nTids; i++){
					final int tid=data.getInt();
					final long first=data.getLong();
					final int count=data.getInt();
					if(tid<1 || count<1 || first!=priorEnd || first+count>genes ||
							tidToEntry.contains(tid)){
						throw new IllegalArgumentException("Malformed contribution directory entry "+i+".");
					}
					tids[i]=tid; firstGenes[i]=first; counts[i]=count;
					tidToEntry.put(tid,i); priorEnd=first+count;
				}
				if(priorEnd!=genes){throw new IllegalArgumentException("Contribution directory does not cover all genes.");}
				success=true;
				return new Reader(channel,data,genes,tids,firstGenes,counts,tidToEntry);
			}finally{
				if(!success){channel.close();}
			}
		}

		/** Returns the reference-gene count for one TID, or throws if absent. */
		public int nativeTotal(final int tid){return counts[entry(tid)];}

		/** Returns one gene's physical CDS length. */
		public int length(final int tid, final int ordinal){
			final long offset=recordOffset(tid,ordinal);
			return data.getInt((int)offset);
		}

		/** Returns one gene's tracked-family rank, or {@code -1} when unmapped/masked/bait/rejected. */
		public int familyRank(final int tid, final int ordinal){
			final long offset=recordOffset(tid,ordinal)+4;
			return data.getChar((int)offset)-1;
		}

		/** Number of indexed TIDs. */
		public int tids(){return tids.length;}

		/** Number of indexed gene identities. */
		public long genes(){return genes;}

		private int entry(final int tid){
			final int entry=tidToEntry.get(tid);
			if(entry<0){throw new IllegalArgumentException("TID absent from contribution index: "+tid);}
			return entry;
		}

		private long recordOffset(final int tid, final int ordinal){
			final int entry=entry(tid),count=counts[entry];
			if(ordinal<0 || ordinal>=count){
				throw new IllegalArgumentException("Gene ordinal "+ordinal+" outside [0,"+
					count+") for tid "+tid+".");
			}
			return HEADER_BYTES+(firstGenes[entry]+ordinal)*RECORD_BYTES;
		}

		@Override
		public void close() throws IOException{channel.close();}
	}

	private static void build(final Config c) throws Exception{
		verifySha256(c.manifest,c.manifestSha256,"contribution manifest");
		verifySha256(c.cache,c.cacheSha256,"whole-genome cache");
		final ArrayList<Table> tables=loadManifest(c);
		final File output=new File(c.index).getAbsoluteFile();
		if(output.exists() && !c.overwrite){throw new IllegalArgumentException("index= exists: "+output);}
		final File parent=output.getParentFile();
		if(parent!=null && !parent.isDirectory()){
			throw new IllegalArgumentException("Index parent directory does not exist: "+parent);
		}
		final File partial=File.createTempFile(output.getName()+".",".partial",parent);
		boolean published=false;
		try{
			final BuildState state=writeIndex(partial,tables,c);
			verifySha256(c.manifest,c.manifestSha256,"contribution manifest changed during build");
			verifySha256(c.cache,c.cacheSha256,"whole-genome cache changed during build");
			ReferenceCdsShredSurvivalTable.publishCompleted(partial,output,c.overwrite);
			published=true;
			System.err.println("ReferenceCdsGeneContributionIndex PASS: tables="+
				tables.size()+" tids="+state.entries.size()+" genes="+state.genes+
				" index="+output);
		}finally{if(!published && partial.exists()){partial.delete();}}
	}

	private static BuildState writeIndex(final File partial,
			final ArrayList<Table> tables, final Config c) throws Exception{
		final BuildState state=new BuildState();
		for(final Table table : tables){scanTable(table,c,state,null,true);}
		if(state.genes!=c.expectedGenes){
			throw new IllegalArgumentException("Indexed genes "+state.genes+
				" != expected "+c.expectedGenes+".");
		}
		long first=0;
		for(final Entry entry : state.entries){
			if(entry.count<1){throw new AssertionError("Contribution directory contains an empty TID.");}
			entry.firstGene=first; first+=entry.count;
		}
		if(first!=state.genes){throw new AssertionError("Directory gene total mismatch.");}
		final long directoryOffset=HEADER_BYTES+state.genes*RECORD_BYTES;
		final long fileBytes=directoryOffset+(long)state.entries.size()*DIRECTORY_BYTES;
		if(fileBytes>Integer.MAX_VALUE){
			throw new IllegalArgumentException("Contribution index exceeds schema-1 mmap range: "+fileBytes);
		}
		try(FileChannel channel=FileChannel.open(partial.toPath(),StandardOpenOption.READ,
				StandardOpenOption.WRITE,StandardOpenOption.TRUNCATE_EXISTING)){
			channel.position(fileBytes-1);
			writeFully(channel,ByteBuffer.wrap(new byte[1]));
			final MappedByteBuffer data=channel.map(FileChannel.MapMode.READ_WRITE,0,fileBytes);
			data.order(ByteOrder.BIG_ENDIAN);
			data.put(MAGIC).putInt(VERSION).putInt(HEADER_BYTES).putInt(c.expectedFamilies)
				.putInt(state.entries.size()).putLong(state.genes).putLong(directoryOffset)
				.put(unhex(c.manifestSha256)).put(unhex(c.cacheSha256));
			state.writeOrdinals=new int[state.entries.size()];
			for(final Table table : tables){scanTable(table,c,state,data,false);}
			for(int i=0; i<state.entries.size(); i++){
				if(state.writeOrdinals[i]!=state.entries.get(i).count){
					throw new AssertionError("Second-pass contribution count mismatch for tid "+
						state.entries.get(i).tid+".");
				}
			}
			int position=(int)directoryOffset;
			for(final Entry entry : state.entries){
				data.putInt(position,entry.tid); position+=4;
				data.putLong(position,entry.firstGene); position+=8;
				data.putInt(position,entry.count); position+=4;
			}
			if(position!=fileBytes){throw new AssertionError("Contribution directory byte-count mismatch.");}
			data.force(); channel.force(true);
		}
		return state;
	}

	private static void scanTable(final Table table, final Config c,
			final BuildState state, final MappedByteBuffer data,
			final boolean countPass) throws Exception{
		final File source=new File(table.path);
		final long sourceBytes=source.length(),sourceModified=source.lastModified();
		if(sourceBytes<1 || sourceModified<1){
			throw new IllegalArgumentException("Unreadable contribution table: "+table.path);
		}
		final ByteFile input=ByteFile.makeByteFile(table.path,false);
		final LineParser1 parser=new LineParser1((byte)'\t');
		boolean schema=false,columns=false,end=false;
		long rows=0,mapped=0,masked=0,bait=0;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				if(line.length==0){throw new IllegalArgumentException("Blank contribution row: "+table.path);}
				parser.set(line);
				if(line[0]=='#'){
					if(parser.termEquals("#schema_version",0)){
						if(schema || parser.terms()!=2 || !parser.termEquals(ReferenceCdsGeneContributionTable.SCHEMA,1)){
							throw new IllegalArgumentException("Bad contribution schema: "+table.path);
						}
						schema=true;
					}else if(parser.termEquals("#expected_families",0)){
						if(parser.terms()!=2 || parser.parseInt(1)!=c.expectedFamilies){
							throw new IllegalArgumentException("Contribution family count mismatch: "+table.path);
						}
					}else if(parser.termEquals("#columns",0)){
						final String expected="#columns\t"+ReferenceCdsGeneContributionTable.COLUMNS;
						if(columns || !expected.equals(new String(line,StandardCharsets.UTF_8))){
							throw new IllegalArgumentException("Contribution columns mismatch: "+table.path);
						}
						columns=true;
					}else if(parser.termEquals("#end",0)){
						end=true;
					}
					continue;
				}
				if(!schema || !columns || end || parser.terms()!=11 || !parser.termEquals("G",0)){
					throw new IllegalArgumentException("Malformed contribution data row in "+table.path);
				}
				final int tid=parser.parseInt(2),ordinal=parser.parseInt(3);
				final int length=parser.parseInt(8);
				if(tid<1 || length<1){throw new IllegalArgumentException("Invalid contribution tid/length in "+table.path);}
				int entryIndex=state.tidToEntry.get(tid);
				if(countPass && entryIndex<0){
					entryIndex=state.entries.size();
					state.tidToEntry.put(tid,entryIndex);
					state.entries.add(new Entry(tid));
				}
				if(entryIndex<0){throw new AssertionError("Second-pass TID absent from first pass: "+tid);}
				final Entry entry=state.entries.get(entryIndex);
				final int expectedOrdinal=countPass ? entry.count : state.writeOrdinals[entryIndex];
				if(ordinal!=expectedOrdinal){
					throw new IllegalArgumentException("Schema 1 requires one physical row per contiguous gene ordinal; tid="+
						tid+" expected="+expectedOrdinal+" observed="+ordinal+".");
				}
				final int family;
				if(parser.termEquals("ASSIGNED",10)){
					family=parser.parseInt(9);
					if(family<0 || family>=c.expectedFamilies){throw new IllegalArgumentException("Assigned family outside roster.");}
					mapped++;
				}else{
					if(!parser.termEquals("-",9) || !knownNonassignedOutcome(parser)){
						throw new IllegalArgumentException("Invalid non-assigned contribution outcome in "+table.path);
					}
					family=-1;
					if(parser.termEquals("MASKED_INTERCEPTED",10)){masked++;}
					else if(parser.termEquals("BAIT_INTERCEPTED",10)){bait++;}
				}
				if(countPass){
					entry.count++; state.genes++;
				}else{
					final long record=entry.firstGene+ordinal;
					final int offset=(int)(HEADER_BYTES+record*RECORD_BYTES);
					data.putInt(offset,length); data.putChar(offset+4,(char)(family+1));
					state.writeOrdinals[entryIndex]++;
				}
				rows++;
			}
		}finally{input.close();}
		if(!schema || !columns || !end || rows!=table.rows || rows!=table.identities ||
				mapped!=table.mapped || masked!=table.masked || bait!=table.bait){
			throw new IllegalArgumentException("Contribution manifest/table count mismatch for "+table.phylum+".");
		}
		if(countPass){verifySha256(table.path,table.sha256,"contribution table "+table.phylum);}
		if(source.length()!=sourceBytes || source.lastModified()!=sourceModified){
			throw new IllegalArgumentException("Contribution table changed during indexing: "+table.path);
		}
	}

	private static boolean knownNonassignedOutcome(final LineParser1 parser){
		return parser.termEquals("MASKED_INTERCEPTED",10) ||
			parser.termEquals("BAIT_INTERCEPTED",10) ||
			parser.termEquals("NO_PRESELECT_SURVIVOR",10) ||
			parser.termEquals("NO_ACCEPTANCE_PASS",10) ||
			parser.termEquals("NO_VALID_KMER",10) ||
			parser.termEquals("INVALID_SEQUENCE",10);
	}

	private static ArrayList<Table> loadManifest(final Config c){
		final ArrayList<Table> tables=new ArrayList<Table>();
		final ByteFile input=ByteFile.makeByteFile(c.manifest,false);
		final LineParser1 parser=new LineParser1((byte)'\t');
		boolean schema=false,columns=false,end=false;
		long rows=0,identities=0,mapped=0,masked=0,bait=0;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				if(line.length==0){throw new IllegalArgumentException("Blank contribution-manifest row.");}
				parser.set(line);
				if(line[0]=='#'){
					if(parser.termEquals("#schema_version",0)){
						if(schema || parser.terms()!=2 ||
								!parser.termEquals("reference_cds_gene_contribution_manifest_v1",1)){
							throw new IllegalArgumentException("Contribution-manifest schema mismatch.");
						}
						schema=true;
					}else if(parser.termEquals("#columns",0)){
						if(columns || !MANIFEST_COLUMNS.equals(new String(line,StandardCharsets.UTF_8))){
							throw new IllegalArgumentException("Contribution-manifest columns mismatch.");
						}
						columns=true;
					}else if(parser.termEquals("#end",0)){end=true;}
					continue;
				}
				if(!schema || !columns || end || parser.terms()!=9){
					throw new IllegalArgumentException("Malformed contribution-manifest data row.");
				}
				final String phylum=parser.parseString(0),path=parser.parseString(2),sha256=parser.parseString(3);
				final long tableRows=parser.parseLong(4),tableIdentities=parser.parseLong(5),
					tableMapped=parser.parseLong(6),tableMasked=parser.parseLong(7),tableBait=parser.parseLong(8);
				if(phylum.length()==0 || path.length()==0 || !isHex64(sha256) || tableRows<1 ||
						tableIdentities!=tableRows || tableMapped<0 || tableMasked<0 || tableBait<0 ||
						tableMapped+tableMasked+tableBait>tableRows){
					throw new IllegalArgumentException("Invalid contribution-manifest row for "+phylum+".");
				}
				tables.add(new Table(phylum,path,sha256,tableRows,tableIdentities,
					tableMapped,tableMasked,tableBait));
				rows+=tableRows; identities+=tableIdentities; mapped+=tableMapped;
				masked+=tableMasked; bait+=tableBait;
			}
		}finally{input.close();}
		if(!schema || !columns || !end || tables.size()!=c.expectedTables ||
				rows!=c.expectedGenes || identities!=rows || mapped!=c.expectedMapped ||
				masked!=c.expectedMasked || bait!=c.expectedBait){
			throw new IllegalArgumentException("Contribution-manifest global totals mismatch.");
		}
		return tables;
	}

	private static void writeFully(final FileChannel channel,
			final ByteBuffer buffer) throws IOException{
		while(buffer.hasRemaining()){channel.write(buffer);}
	}

	private static void verifySha256(final String path, final String expected,
			final String label){
		requireHex64(expected,label+" expected SHA-256");
		final String observed=ReferenceCdsSurvivalLabelReader.sha256File(path);
		if(!expected.equals(observed)){
			throw new IllegalArgumentException(label+" SHA-256 mismatch: expected_sha80="+
				expected.substring(44)+" observed_sha80="+observed.substring(44)+" for "+path);
		}
	}

	private static void requireHex64(final String text, final String label){
		if(!isHex64(text)){throw new IllegalArgumentException(label+" must be 64 lowercase hex characters.");}
	}

	private static boolean isHex64(final String text){
		if(text==null || text.length()!=64){return false;}
		for(int i=0; i<64; i++){
			final char c=text.charAt(i);
			if(!((c>='0' && c<='9') || (c>='a' && c<='f'))){return false;}
		}
		return true;
	}

	private static byte[] unhex(final String text){
		requireHex64(text,"SHA-256");
		final byte[] bytes=new byte[32];
		for(int i=0; i<32; i++){
			bytes[i]=(byte)((Character.digit(text.charAt(2*i),16)<<4) |
				Character.digit(text.charAt(2*i+1),16));
		}
		return bytes;
	}

	private static String hex(final byte[] bytes){
		final char[] chars=new char[bytes.length*2];
		for(int i=0; i<bytes.length; i++){
			final int value=bytes[i]&0xff;
			chars[2*i]=HEX[value>>>4]; chars[2*i+1]=HEX[value&15];
		}
		return new String(chars);
	}

	private static final class BuildState{
		long genes=0;
		final ArrayList<Entry> entries=new ArrayList<Entry>(1<<16);
		final IntHashMap tidToEntry=new IntHashMap(1<<16);
		int[] writeOrdinals;
	}

	private static final class Entry{
		final int tid;
		int count;
		long firstGene;
		Entry(final int tid_){tid=tid_;}
	}

	private static final class Table{
		final String phylum,path,sha256;
		final long rows,identities,mapped,masked,bait;
		Table(final String phylum_, final String path_, final String sha256_,
				final long rows_, final long identities_, final long mapped_,
				final long masked_, final long bait_){
			phylum=phylum_; path=path_; sha256=sha256_; rows=rows_;
			identities=identities_; mapped=mapped_; masked=masked_; bait=bait_;
		}
	}

	private static final class Config{
		String manifest,manifestSha256,cache,cacheSha256,index,indexSha256;
		String survivalManifest,survivalSha80,shredsCache,shredsSha80;
		int expectedFamilies,expectedTables;
		long expectedGenes,expectedMapped,expectedMasked,expectedBait;
		boolean build,overwrite;

		static Config parse(final String[] args){
			final HashMap<String,String> map=new HashMap<String,String>();
			for(final String arg : args){
				final int equals=arg.indexOf('=');
				if(equals<1){throw new IllegalArgumentException("Expected key=value: "+arg);}
				final String key=arg.substring(0,equals).toLowerCase(Locale.ROOT);
				if(map.put(key,arg.substring(equals+1))!=null){throw new IllegalArgumentException("Duplicate argument: "+key);}
			}
			final Config c=new Config();
			// Alias exclusion precedes file reads, even when a build does not use an index pin.
			for(final String stem : new String[]{"manifest", "cache", "index"}){
				if(map.containsKey(stem+"sha256") && map.containsKey(stem+"sha80")){
					throw new IllegalArgumentException("Specify only one of "+stem+"sha256= or "+stem+"sha80=");
				}
			}
			c.index=req(map,"index"); c.indexSha256=map.get("indexsha256");
			c.manifest=req(map,"manifest"); c.manifestSha256=inputSha(map, "manifest", c.manifest);
			c.cache=req(map,"cache"); c.cacheSha256=inputSha(map, "cache", c.cache);
			c.expectedFamilies=positiveInt(map,"expectedfamilies");
			if(c.expectedFamilies>Character.MAX_VALUE){
				throw new IllegalArgumentException("Schema 1 stores family rank+1 as an unsigned short; expectedfamilies must be <= "+
					(int)Character.MAX_VALUE+": "+c.expectedFamilies);
			}
			c.build=parse.Parse.parseBoolean(optional(map,"build","false"));
			if(map.containsKey("survivalmanifest") || map.containsKey("survivalsha80") ||
				map.containsKey("shredscache") || map.containsKey("shredssha80")){
				if(c.build){throw new IllegalArgumentException("Survival count comparison requires build=f");}
				c.survivalManifest=req(map, "survivalmanifest");
				c.survivalSha80=req(map, "survivalsha80");
				c.shredsCache=req(map, "shredscache");
				c.shredsSha80=req(map, "shredssha80");
				if(!c.survivalSha80.matches("[0-9a-f]{20}") || !c.shredsSha80.matches("[0-9a-f]{20}")){
					throw new IllegalArgumentException("survivalsha80 and shredssha80 require 20 lowercase hexadecimal characters");
				}
			}
			if(map.containsKey("overwrite") && map.containsKey("ow")){
				throw new IllegalArgumentException("Specify only one of overwrite= or ow=.");
			}
			c.overwrite=parse.Parse.parseBoolean(map.containsKey("ow") ? map.get("ow") :
				optional(map,"overwrite","false"));
			if(c.build){
				if(map.containsKey("indexsha80")){throw new IllegalArgumentException("indexsha80 requires build=f; index is an output during build=t");}
				c.expectedTables=positiveInt(map,"expectedtables");
				c.expectedGenes=positiveLong(map,"expectedgenes");
				c.expectedMapped=nonnegativeLong(map,"expectedmapped");
				c.expectedMasked=nonnegativeLong(map,"expectedmasked");
				c.expectedBait=nonnegativeLong(map,"expectedbait");
			}else{
				c.indexSha256=inputSha(map, "index", c.index);
			}
			for(final String key : map.keySet()){
				if(!ALLOWED.containsKey(key)){throw new IllegalArgumentException("Unknown argument: "+key);}
			}
			return c;
		}

		private static String req(final HashMap<String,String> map, final String key){
			final String value=map.get(key);
			if(value==null || value.length()==0){throw new IllegalArgumentException("Required: "+key+"=");}
			return value;
		}
		private static String optional(final HashMap<String,String> map,
				final String key, final String defaultValue){
			final String value=map.get(key); return value==null ? defaultValue : value;
		}
		/** Resolves one operator pin to the full digest required by legacy loaders and index metadata. */
		private static String inputSha(final HashMap<String,String> map, final String stem, final String path){
			assert(path!=null && !path.isEmpty()) : "Config.req validates the input path before digest resolution";
			final String shortKey=stem+"sha80", longKey=stem+"sha256";
			if(map.containsKey(shortKey)){
				final String pin=req(map, shortKey);
				if(!pin.matches("[0-9a-f]{20}")){throw new IllegalArgumentException(shortKey+" must be 20 lowercase hex characters");}
				return verifySha80(path, pin);
			}
			return requiredSha(req(map, longKey), longKey);
		}
		private static String requiredSha(final String value, final String key){
			requireHex64(value,key); return value;
		}
		private static int positiveInt(final HashMap<String,String> map, final String key){
			final long value=positiveLong(map,key);
			if(value>Integer.MAX_VALUE){throw new IllegalArgumentException(key+" overflow.");}
			return (int)value;
		}
		private static long positiveLong(final HashMap<String,String> map, final String key){
			final long value=nonnegativeLong(map,key);
			if(value<1){throw new IllegalArgumentException(key+" must be positive.");}
			return value;
		}
		private static long nonnegativeLong(final HashMap<String,String> map,
				final String key){
			final String text=req(map,key);
			try{
				final long value=Long.parseLong(text);
				if(value<0 || !Long.toString(value).equals(text)){throw new NumberFormatException();}
				return value;
			}catch(NumberFormatException e){
				throw new IllegalArgumentException(key+" must be a canonical nonnegative long.",e);
			}
		}
	}

	private static final HashMap<String,Boolean> ALLOWED=new HashMap<String,Boolean>();
	static{
		for(final String key : new String[]{"manifest","manifestsha256","cache",
			"cachesha256","index","indexsha256","expectedfamilies","expectedtables",
			"expectedgenes","expectedmapped","expectedmasked","expectedbait","build",
			"overwrite","ow","survivalmanifest","survivalsha80","shredscache","shredssha80",
			"manifestsha80","cachesha80","indexsha80"}){ALLOWED.put(key,Boolean.TRUE);}
	}

	private static final byte[] MAGIC={'M','A','G','Q','C','G','I','1'};
	private static final int VERSION=1,HEADER_BYTES=104;
	private static final int RECORD_BYTES=6,DIRECTORY_BYTES=16;
	private static final char[] HEX="0123456789abcdef".toCharArray();
	private static final String MANIFEST_COLUMNS="#columns\tgff_phylum\tassignment_phylum\t"+
		"table_path\ttable_sha256\trows\tidentities\tmapped\tmasked\tbait";
}
