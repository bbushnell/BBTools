package dna;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.concurrent.ConcurrentHashMap;

import align2.AbstractIndex;
import align2.AbstractMapper;
import align2.BBSplitter;
import align2.ChromLoadThread;
import align2.RefToIndex;
import driver.Search;
import fileIO.ChainBlock;
import fileIO.ChainLine;
import fileIO.ReadWrite;
import fileIO.TextFile;
import parse.Parse;
import server.PercentEncoding;
import shared.Primes;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.Range;

/**
 * Process-wide reference sequences, scaffold metadata, legacy gene annotations,
 * resource paths and external-tool availability caches.
 * Genome and path initialization must finish before publishing data to workers;
 * callers must coordinate resets and mutation of the returned shared arrays/maps.
 * Loader locks serialize construction but do not make every unsynchronized cache
 * read or public field mutation safe. Scaffold-name publication has its own
 * volatile guarantee, described in {@link #scafNameTable()}.
 * Chromosome positions and range endpoints are zero-based and inclusive.
 *
 * @author Brian Bushnell
 * @date 2014
 */
public class Data{

	/** Empty main method for class testing purposes.
	 * @param args Command-line arguments (unused) */
	public static void main(String[] args){}

	//TODO IMPORTANT!  Ensure that this unloads everything big, AND that reloading subsequently works OK.
	/**
	 * Clears chromosome/index/scaffold and active-build gene caches, mapping bounds and BBSplitter tables.
	 * Fixed gene_names_37.1 dictionaries, chain maps, paths, interned strings and tool probes remain.
	 * Outer gene-cache arrays retain their identities; previously borrowed inner objects are not mutated.
	 * This is not a complete initial-state reset and must not race with readers.
	 */
	public static void unloadAll(){
		chromosomePlusMatrix=null;
		AbstractIndex.clear();
		RefToIndex.clear();

		AbstractMapper.minChrom=1;
		AbstractMapper.maxChrom=Integer.MAX_VALUE;

		numChroms=0;
		numBases=0;
		numDefinedBases=0;
		numContigs=0;
		numScaffolds=0;
		interScaffoldPadding=0;
		chromLengths=null;
		chromDefinedBases=null;
		chromUndefinedBases=null;
		chromContigs=null;
		chromScaffolds=null;
		chromStartPad=null;

		scaffoldNames=null;
		scaffoldLocs=null;
		scaffoldLengths=null;

		BBSplitter.setCountTable=null;
		BBSplitter.scafCountTable=null;
		BBSplitter.streamTable=null;
		BBSplitter.streamTableAmbiguous=null;

		scaffoldNameTable=null;
		genomeSource=null;
		name=null;

		//DNA-027: clear every chromosome slot so independent lazy getters cannot return the previous build.
		Arrays.fill(geneMatrix, null);
		Arrays.fill(geneSetMatrix, null);
		Arrays.fill(geneTxRangeMatrix, null);
		Arrays.fill(geneSetRangeMatrix, null);
		Arrays.fill(geneCodeRangeMatrix, null);
		Arrays.fill(geneCodeAndExonRangeMatrix, null);
		Arrays.fill(geneNearbyRangeMatrix, null);
		Arrays.fill(exonRangeMatrix, null);
		Arrays.fill(geneNameTable, null);
		Arrays.fill(transcriptNameTable, null);
		geneIDTable=null;

		GENOME_BUILD=-1;
		genome_set_to=-1;
	}

	//TODO IMPORTANT!  Ensure that this unloads everything big, AND that reloading subsequently works OK.
	/**
	 * Drops only the cached sequence for a chromosome; existing borrowed references remain valid.
	 * @param chrom Chromosome number to unload
	 * @param unloadSoft Whether to perform soft unload (currently unused)
	 */
	public static void unload(int chrom, boolean unloadSoft){

//		unloadGenes(chrom);

		chromosomePlusMatrix[chrom]=null;
	}

	/**
	 * Drops per-chromosome gene, gene-set and range arrays.
	 * Name/transcript lookup tables and the active gene-ID map are retained (DNA-053).
	 * Callers must exclude concurrent readers and loaders.
	 * @param chrom Chromosome number to unload gene data for
	 */
	public static void unloadGenes(int chrom){
		//TODO: Probable bug - DNA-053: retained symbol/transcript/active-ID lookups may still return this chromosome's old annotations.
		geneMatrix[chrom]=null;
		geneSetMatrix[chrom]=null;
		geneTxRangeMatrix[chrom]=null;
		geneSetRangeMatrix[chrom]=null;
		geneCodeRangeMatrix[chrom]=null;
		geneCodeAndExonRangeMatrix[chrom]=null;
		geneNearbyRangeMatrix[chrom]=null;
		exonRangeMatrix[chrom]=null;
	}

	/**
	 * Finds the index of a value in a byte array.
	 * @param x Value to search for
	 * @param array Array to search in
	 * @return First index narrowed to byte, or -1 if absent; see DNA-025 for overflow
	 */
	public static byte find(int x, byte[] array){
		//TODO: Probable bug - DNA-025: indices >=128 wrap in the byte return; index255 collides with the absent sentinel. No live caller confirmed.
		for(int i=0; i<array.length; i++){ //int counter: byte i overflowed (AIOOBE) for arrays >=128
			if(array[i]==x){return (byte)i;}
		}
		return -1;
	}

	/** Reverses a byte array in-place.
	 * @param array Array to reverse */
	public static void reverse(byte[] array){
		int mid=array.length/2;
		for(int i=0; i<mid; i++){
			byte temp=array[i];
			array[i]=array[array.length-i-1];
			array[array.length-i-1]=temp;
		}
	}

	/**
	 * Retrieves all genes for specified chromosome, loading if necessary.
	 * Returns the shared cached array; see the class initialization/threading contract.
	 * @param chrom Chromosome number
	 * @return Array of genes on the chromosome
	 */
	public static Gene[] getGenes(int chrom){
		if(geneMatrix[chrom]==null){
			loadGenes(chrom);
		}
		return geneMatrix[chrom];
	}

	/**
	 * Retrieves genes for specified chromosome filtered by strand.
	 * @param chrom Chromosome number
	 * @param strand Strand orientation (+ or -)
	 * @return New array containing shared Gene objects on the specified strand
	 */
	public static Gene[] getGenes(int chrom, byte strand){
		ArrayList<Gene> genes=new ArrayList<Gene>();
		for(Gene g : getGenes(chrom)){
			if(g.strand==strand){
				genes.add(g);
			}
		}
		return genes.toArray(new Gene[genes.size()]);
	}

	/**
	 * Retrieves gene sets for specified chromosome, loading if necessary.
	 * Gene sets group related transcript variants under common identifiers.
	 * @param chrom Chromosome number
	 * @return Shared cached array of gene sets on the chromosome
	 */
	public static GeneSet[] getGeneSets(int chrom){
		if(geneSetMatrix[chrom]==null){
			loadGenes(chrom);
		}
		return geneSetMatrix[chrom];
	}

	/**
	 * Builds and caches a lookup table mapping gene IDs to gene sets.
	 * Scans legacy chromosomes 1-25, independent of numChroms. The lock serializes
	 * construction; the nonvolatile fast-path read requires caller-managed publication.
	 * @return Shared mutable map of gene IDs to lists of shared gene sets
	 */
	public static HashMap<Integer, ArrayList<GeneSet>> getGeneIDTable(){
		if(geneIDTable==null){

//			System.err.println("WAITING FOR CS");
			synchronized(GENEIDLOCK){
//				System.err.println("ENTER CS");
				if(geneIDTable==null){
//					System.err.println("ENTER CS2");
					HashMap<Integer, ArrayList<GeneSet>> temp=new HashMap<Integer, ArrayList<GeneSet>>(2048);
					for(byte chrom=1; chrom<=25; chrom++){
						GeneSet[] set=getGeneSets(chrom);
						for(GeneSet gs : set){
							int id=-1;
							for(Gene g : gs.genes){
								if(id==-1){id=g.id;}
								else{assert(id==g.id) : gs+"\n"+gs.genes+"\n";}
							}
							assert(id>-1);
							Integer key=Integer.valueOf(id);
							ArrayList<GeneSet> value=temp.get(key);
							//					assert(old==null || chrom>22) : "\nCollision!\n\n"+gs+"\n\nis overwriting\n\n"+old;
							if(value==null){
								value=new ArrayList<GeneSet>(2);
								temp.put(key, value);
							}
							value.add(gs);
						}
					}
//					System.err.println("EXIT CS2");
					geneIDTable=temp;
				}
//				System.err.println("EXIT CS");
			}

		}
//		System.err.println("GeneIDTable contains "+geneIDTable.size()+" entries.");
		return geneIDTable;
	}

	/**
	 * Retrieves chromosome sequence data, loading from disk if necessary.
	 * Uses striped loader locks; see the class publication/reset contract.
	 * Asserts that base zero is N, a required leading-padding check.
	 *
	 * @param chrom Chromosome number
	 * @return Shared ChromosomeArray containing sequence data
	 */
	public static ChromosomeArray getChromosome(int chrom){
		assert(chromosomePlusMatrix!=null);
		assert(chromosomePlusMatrix.length>chrom) : chromosomePlusMatrix.length+", "+chrom;
		if(chromosomePlusMatrix[chrom]==null){
			synchronized(CHROMLOCKS[chrom%CHROMLOCKS.length]){
				if(chromosomePlusMatrix[chrom]==null){loadChromosome(chrom);}
			}
		}
		assert(chromosomePlusMatrix[chrom].array[0]=='N') : (char)chromosomePlusMatrix[chrom].array[0]+
			"\nIf you see this message, please regenerate your index.\n"/*+new String(chromosomePlusMatrix[chrom].array)*/;//startpad was too low or for some reason invalid.
		return chromosomePlusMatrix[chrom];
	}

	/**
	 * Loads gene data for specified chromosome from serialized files.
	 * Creates gene matrices, range matrices, and name lookup tables.
	 * Handles both transcript-based and gene symbol-based indexing.
	 * @param chrom Chromosome number to load
	 */
	private static void loadGenes(int chrom){

		if(geneMatrix[chrom]!=null){return;} //In case another thread already loaded the chromosome
		synchronized(CHROMLOCKS[chrom%CHROMLOCKS.length]){
			if(geneMatrix[chrom]==null){

				//		Gene[] genes=FindExons.readGenes(ROOT_GENE+"ref/chr"+chrom+".Ref.Table", Gene.FORMAT_NM);
				//		Gene[] genes=FindExons.readGenes(ROOT_GENE+"ref2/ccds-chr"+chrom+"-genes.txt", Gene.FORMAT_CCDS);
				//		Gene[] genes=FindExons.readGenes(ROOT_GENE+"ref3/ccds-chr"+chrom+"-genes.txt", Gene.FORMAT_CCDS);

				//		Gene[] genes=ReadWrite.readArray(Gene.class, ROOT_GENE+"seqGene/chr"+chrom+".ga");

				String fname=ROOT_GENE+"Build"+GENOME_BUILD+"/"+GENE_MAP+"/chr"+chrom+".ga";

				Gene[] genes=ReadWrite.readArray(Gene.class, fname, true);

				Arrays.sort(genes);
//				geneMatrix[chrom]=genes;

				geneTxRangeMatrix[chrom]=findGeneRanges(genes, TX_RANGE);
				//TODO: Probable bug - DNA-036: compareTo sorts by txStart, not codeStart; CODE_RANGE can emit unsorted/unmerged ranges.
				//For example, txStart/codeStart pairs0/100 then10/20 reverse coding order. Establish live consumers before repair.
				geneCodeRangeMatrix[chrom]=findGeneRanges(genes, CODE_RANGE);
				geneCodeAndExonRangeMatrix[chrom]=findCodeAndExonRanges(genes, false, true);
				exonRangeMatrix[chrom]=findCodeAndExonRanges(genes, false, false);
				geneNearbyRangeMatrix[chrom]=findCodeAndExonRanges(genes, true, true);

				HashMap<String, ArrayList<Gene>> temp=new HashMap<String, ArrayList<Gene>>();
				HashMap<String, GeneSet> gntable=new HashMap<String, GeneSet>();
				HashMap<String, Gene> tntable=new HashMap<String, Gene>();

				for(Gene g : genes){

					String trkey=g.mrnaAcc;
					if(trkey==null){trkey=g.chromosome+"_"+g.id;}
					if(trkey!=null){
						Gene old=tntable.get(trkey);
						if(old!=null){
							//					stdout.println("For transcript '"+g.nameTranscript+"': Overwrote \n"+old+"\nwith\n"+g+"\n");
						}
						tntable.put(trkey, g);
					}

					String key=g.symbol;
					if(key==null){key=g.mrnaAcc;}
					ArrayList<Gene> list=temp.get(key);
					if(list==null){
						list=new ArrayList<Gene>(8);
						temp.put(key, list);
					}
					list.add(g);
				}

				GeneSet[] gsm=new GeneSet[temp.size()];
				String[] keys=temp.keySet().toArray(new String[temp.size()]);
				for(int i=0; i<keys.length; i++){
					String key=keys[i];
					ArrayList<Gene> list=temp.get(key);
					GeneSet gs=new GeneSet(key, list);
					gsm[i]=gs;
					gntable.put(key, gs);
				}

				geneNameTable[chrom]=gntable;
				transcriptNameTable[chrom]=tntable;
				geneSetMatrix[chrom]=gsm;
				Arrays.sort(geneSetMatrix[chrom]);

				geneSetRangeMatrix[chrom]=findGeneSetRanges(geneSetMatrix[chrom]);

				{
					assert(geneMatrix[chrom]==null) : "Need to sync.";
					geneMatrix[chrom]=genes;
				}
			}
		}
	}

	/**
	 * Loads chromosome sequence data for range of chromosomes.
	 * Uses ChromLoadThread for parallel loading into the existing cache.
	 * The bulk-load monitor differs from the striped lazy-loader monitors;
	 * callers must coordinate bulk loading with lazy access and resets.
	 * @param min Minimum chromosome number to load
	 * @param max Maximum chromosome number to load
	 */
	public static void loadChromosomes(int min, int max){
		synchronized(CHROMLOCKS){
			String pattern=chromFname(GENOME_BUILD);
			ChromLoadThread.loadAll(pattern, min, max, chromosomePlusMatrix);
		}
	}

	/**
	 * Loads single chromosome sequence data from disk.
	 * Reads serialized ChromosomeArray and validates chromosome number.
	 * @param chrom Chromosome number to load
	 */
	private static void loadChromosome(int chrom){
//		assert(false);
		assert(chromosomePlusMatrix[chrom]==null);
//		assert(chrom>0) : chrom; //No longer valid since chrom 0 is now semi-allowed

		String fname=chromFname(chrom, GENOME_BUILD);
		sysout.println("Loading "+fname);
		chromosomePlusMatrix[chrom]=ReadWrite.read(ChromosomeArray.class, fname, false);
		assert(chromosomePlusMatrix[chrom].chromosome==chrom);
	}

	/** Returns file extension for chromosome files (.chrom or .chrom.gz) */
	public static final String chromExtension(){
		return ".chrom"+(CHROMGZ ? ".gz" : "");
	}

	/**
	 * Constructs file path for specific chromosome in genome build.
	 * @param chrom Chromosome number
	 * @param genome Genome build number
	 * @return Full file path to chromosome data file
	 */
	public static final String chromFname(int chrom, int genome){
		return ROOT_GENOME+genome+"/chr"+chrom+chromExtension();
	}

	/**
	 * Constructs file path pattern for all chromosomes in genome build.
	 * Uses '#' placeholder for chromosome number.
	 * @param genome Genome build number
	 * @return File path pattern with chromosome placeholder
	 */
	public static final String chromFname(int genome){
		return ROOT_GENOME+genome+"/chr#"+chromExtension();
	}

	/**
	 * Creates merged ranges from overlapping or adjacent gene coordinates.
	 * Groups genes into continuous genomic regions based on transcription
	 * or coding boundaries. Consumes the supplied order without sorting; callers
	 * must provide an order suitable for the selected coordinates. Each result's
	 * obj1 contains a new array of shared Gene references.
	 *
	 * @param genes Array of genes to process
	 * @param mode Range type (TX_RANGE for transcription, CODE_RANGE for coding)
	 * @return Array of merged ranges containing grouped genes
	 */
	public static Range[] findGeneRanges(Gene[] genes, final int mode){

		ArrayList<Range> list=new ArrayList<Range>(8192);
		ArrayList<Gene> glist=new ArrayList<Gene>(64);

		Range current=null;
		for(int i=0; i<genes.length; i++){
			Gene g=genes[i];
			Range r;

			int a, b;

			switch(mode){
				case TX_RANGE: {a=g.txStart; b=g.txStop;}
				break;

				case CODE_RANGE: {a=g.codeStart; b=g.codeStop;}
				break;

				default: {throw new RuntimeException();}
			}

			if(b>=a){
				r=new Range(a, b);

				if(current==null){
					current=r;
					glist.add(g);
				}else if(current.touches(r)){
					current=current.merge(r);
					glist.add(g);
				}else{
					current.obj1=glist.toArray(new Gene[glist.size()]);
					glist.clear();
					glist.add(g);
					list.add(current);
					current=r;
				}
			}
		}
		if(current!=null){ //i.e., if there were any genes
			current.obj1=glist.toArray(new Gene[glist.size()]);
			list.add(current);
		}

		return list.toArray(new Range[list.size()]);
	}

	/**
	 * Creates merged ranges from gene set coordinates with nearby padding.
	 * Extends boundaries by NEAR distance to capture nearby regulatory regions.
	 * Consumes sorted input without sorting, merging overlaps and adjacency.
	 * Result obj1 arrays retain shared GeneSet references in input order.
	 * @param genes Array of gene sets to process
	 * @return Array of merged ranges with padded boundaries
	 */
	public static Range[] findGeneSetRanges(GeneSet[] genes){

		ArrayList<Range> list=new ArrayList<Range>(8192);
		ArrayList<GeneSet> glist=new ArrayList<GeneSet>(64);

		Range current=null;
		for(int i=0; i<genes.length; i++){
			GeneSet g=genes[i];
			Range r;

			int a=g.minStart-NEAR, b=g.maxEnd+NEAR;

			if(b>=a){
				r=new Range(a, b);

				if(current==null){
					current=r;
					glist.add(g);
				}else if(current.touches(r)){
					current=current.merge(r);
					glist.add(g);
				}else{
					current.obj1=glist.toArray(new GeneSet[glist.size()]);
					glist.clear();
					glist.add(g);
					list.add(current);
					current=r;
				}
			}
		}
		if(current!=null){ //i.e., if there were any genes
			current.obj1=glist.toArray(new GeneSet[glist.size()]);
			list.add(current);
		}

		return list.toArray(new Range[list.size()]);
	}

	/**
	 * Creates ranges from individual exons with optional coding region filtering.
	 * Can extend boundaries for nearby regions or restrict to coding sequences.
	 * Sorts and merges overlapping or adjacent exon ranges, grouping unique genes
	 * in obj1 arrays with unspecified order. Coding intersection precedes expansion.
	 *
	 * @param genes Array of genes to process exons from
	 * @param nearby Whether to extend boundaries by NEAR distance
	 * @param codingOnly Whether to restrict to coding regions only
	 * @return Array of merged exon ranges with associated genes
	 */
	public static Range[] findCodeAndExonRanges(Gene[] genes, boolean nearby, boolean codingOnly){

		ArrayList<Range> list=new ArrayList<Range>(32768);

		for(int i=0; i<genes.length; i++){
			Gene g=genes[i];
			Range r;

			for(Exon ex : g.exons){

				int a=ex.a, b=ex.b;

				if(codingOnly){
					a=max(ex.a, g.codeStart);
					b=min(ex.b, g.codeStop);
				}

				assert(ex.a<=ex.b);
				assert(g.codeStart<=g.codeStop+1) : g;

				if(nearby){
					a=a-NEAR;
					b=b+NEAR;
				}

				if(a<=b){
					r=new Range(a, b);
					r.obj1=g;
//					r.obj2=ex;
					list.add(r);

				}

			}
		}

		ArrayList<Range> list2=new ArrayList<Range>(list.size());
		Shared.sort(list);

		HashSet<Gene> gset=new HashSet<Gene>(64);
		Range current=null;
		for(Range r : list){
			if(current==null){
				gset.add((Gene)r.obj1);
				current=r;
			}else if(current.touches(r)){
				gset.add((Gene)r.obj1);
				current=current.merge(r);
			}else{
				current.obj1=gset.toArray(new Gene[gset.size()]);
				list2.add(current);
				gset.clear();
				gset.add((Gene)r.obj1);
				current=r;
			}
		}

		if(current!=null){
			current.obj1=gset.toArray(new Gene[gset.size()]);
			list2.add(current);
			Shared.sort(list2);
		}

		return list2.toArray(new Range[list2.size()]);
	}

	/** Returns shared merged gene-set ranges, loading annotations if needed. */
	public static Range[] geneSetRangeMatrix(int chrom){
		if(geneSetRangeMatrix[chrom]==null){
			loadGenes(chrom);
		}
		assert(geneSetRangeMatrix[chrom]!=null);
		return geneSetRangeMatrix[chrom];
	}

	/** Returns shared merged exon ranges, including untranslated exon bases. */
	public static Range[] exonRangeMatrix(int chrom){
		if(exonRangeMatrix[chrom]==null){
			loadGenes(chrom);
		}
		assert(exonRangeMatrix[chrom]!=null);
		return exonRangeMatrix[chrom];
	}

	/** Returns shared ranges formed by intersecting exons with coding coordinates. */
	public static Range[] geneCodeAndExonRangeMatrix(int chrom){
		if(geneCodeAndExonRangeMatrix[chrom]==null){
			loadGenes(chrom);
		}
		assert(geneCodeAndExonRangeMatrix[chrom]!=null);
		return geneCodeAndExonRangeMatrix[chrom];
	}

	/** Returns shared coding-exon ranges expanded by NEAR before merging. */
	public static Range[] geneNearbyRangeMatrix(int chrom){
		if(geneNearbyRangeMatrix[chrom]==null){
			loadGenes(chrom);
		}
		assert(geneNearbyRangeMatrix[chrom]!=null);
		return geneNearbyRangeMatrix[chrom];
	}

	/** Returns the shared symbol-to-gene-set map; missing symbols use transcript accessions. */
	public static HashMap<String, GeneSet> geneNameTable(int chrom){
		if(geneNameTable[chrom]==null){
			loadGenes(chrom);
		}
		assert(geneNameTable[chrom]!=null);
		return geneNameTable[chrom];
	}

	/** Returns the shared accession-to-gene map; missing accessions use chromosome_id. */
	public static HashMap<String, Gene> transcriptNameTable(int chrom){
		if(transcriptNameTable[chrom]==null){
			loadGenes(chrom);
		}
		assert(transcriptNameTable[chrom]!=null);
		return transcriptNameTable[chrom];
	}

	/**
	 * Selects the containing or preceding merged gene-set range by binary search.
	 * Requires a nonempty range array; the selected range need not contain the point.
	 * @param chrom Chromosome number
	 * @param loc Genomic coordinate
	 * @return Shared gene-set array associated with the selected range
	 */
	public static GeneSet[] getNearestGeneSets(int chrom, int loc){
		//TODO: Probable bug - driver/Search#001: an empty annotation range array cannot be indexed by this caller.
		Range[] r=geneSetRangeMatrix(chrom);
		int index=driver.Search.findPointBinary(loc, r);
		GeneSet[] sets=(GeneSet[]) r[index].obj1;
		if(sets==null || sets.length==0){
			assert(false);
			return null;
		}
		return sets;
	}

	/** Selects gene sets for an inclusive interval, or null when no ranges exist.
	 * May return a shared selected-range array or a new array of shared gene sets.
	 * Returns whole payloads of overlapping padded, merged ranges in range order.
	 * A gap selects the closer bounding range by endpoint distance, choosing left on ties;
	 * queries before or after all ranges select the first or last range respectively.
	 * The point overload retains its separate containing-or-preceding behavior. */
	public static GeneSet[] getNearestGeneSets(int chrom, int loc1, int loc2){
		assert(loc2>=loc1);

//		boolean flag=(chrom==21 && loc1<38540895 && loc2>38540895);//TODO UNDO
//
//		if(flag){
//			stdout.println(loc1+", "+loc2+", "+((loc1+loc2)/2));
//			for(GeneSet gs : Data.geneNameTable[chrom].values()){
//				if(gs.intersects(loc1, loc2)){
//					stdout.println("%%% "+gs);
//				}
//			}
//		}

		Range[] ranges=geneSetRangeMatrix(chrom);
		if(ranges==null || ranges.length==0){return null;}
		int index=driver.Search.findPointBinary(loc1, ranges);

//		if(flag){
//			Range r0=ranges[index-1];
//			Range r1=ranges[index];
//			Range r2=ranges[index+1];
//
//			stdout.println("r0: "+r0+"\n"+Arrays.toString((GeneSet[])r0.obj1)+"\n");
//			stdout.println("r1: "+r1+"\n"+Arrays.toString((GeneSet[])r1.obj1)+"\n");
//			stdout.println("r2: "+r2+"\n"+Arrays.toString((GeneSet[])r2.obj1)+"\n");
//
//		}
//		if(flag){stdout.println("c");}

		Range r1=ranges[index];
		Range r2=(index>=ranges.length-1 ? null : ranges[index+1]);

		//DNA-028: Handle gaps before the fast return that also covers nonoverlapping intervals.
		if(r2!=null && loc1>r1.b && loc2<r2.a){
			final long dist1=(long)loc1-r1.b;
			final long dist2=(long)r2.a-loc2;
			if(dist1<=dist2){
				return (GeneSet[])r1.obj1;
			}else{
				return (GeneSet[])r2.obj1;
			}
		}

		if(ranges[index].b>=loc2 || r2==null || r2.a>loc2){
			return (GeneSet[])r1.obj1;
		}

////		if(flag){stdout.println("e");}
//		if(ranges[index].b>=loc2 || (index==ranges.length-1) || ranges[index+1].a>loc2){
////			if(flag){
////				stdout.println("f");
////				stdout.println(ranges[index].b<=loc2);
////				stdout.println((index==ranges.length-1));
////				stdout.println(ranges[index+1].a>loc2);
////				stdout.println(".......");
////			}
//			return sets1;
//		}

//		assert(false) : "Test: This should be very rare, since it is slow.";

		//Otherwise, return all overlapping ranges.
		ArrayList<GeneSet> list=new ArrayList<GeneSet>(4);

		while(index<ranges.length && ranges[index].b<loc1){index++;} //Spin until in range
		for(; index<ranges.length && loc2>=ranges[index].a; index++){
//			if(flag){stdout.println("ADDED RANGE "+ranges[index]);}
			GeneSet[] gsa=(GeneSet[]) ranges[index].obj1;
			for(GeneSet gs : gsa){list.add(gs);}
		}
		return list.toArray(new GeneSet[list.size()]);
	}

	/**
	 * Determines if genomic position falls within exonic regions.
	 *
	 * @param chrom Chromosome number
	 * @param point Genomic coordinate to test
	 * @param thresh Distance threshold for nearby matches
	 * @param isCoding Whether to check coding exons only or all exons
	 * @return true if position is exonic within threshold
	 */
	public static boolean isExonic(byte chrom, int point, int thresh, boolean isCoding){
		//DNA-034: Search owns inclusive threshold-overlap selection; callers must supply a nonnegative thresh.
		Range[] ranges=(isCoding ? Data.geneCodeAndExonRangeMatrix(chrom) : Data.exonRangeMatrix(chrom));
		return Search.containsPointBinary(point, ranges, thresh);
	}

	/**
	 * Prepends whole symbol strings until reaching at least width characters.
	 * Multi-character symbols can overshoot the width; existing text is never truncated.
	 * Returns the original string when no padding is needed, without validating symbol.
	 * A null symbol contributes literal "null" prefixes; a null num still throws.
	 *
	 * @param num String to pad
	 * @param width Target width
	 * @param symbol Prefix string, or null; an empty string is rejected only when padding is required
	 * @return Padded string
	 * @throws IllegalArgumentException if padding is required and symbol is empty
	 * @throws NullPointerException if num is null
	 */
	public static final String padFront(String num, int width, String symbol){
		String r=num;
		//DNA-031: Reject only the no-progress case; preserve all existing finite concatenation behavior.
		if(r.length()<width && symbol!=null && symbol.length()<1){
			throw new IllegalArgumentException("Empty padding symbol would make no progress from length "+r.length()+" to width "+width);
		}
		while(r.length()<width){r=symbol+r;}
		return r;
	}

	/**
	 * Converts number to binary string with specified width.
	 * @param num Number to convert
	 * @param width Minimum width with zero padding
	 * @return Binary string representation
	 */
	public static final String toBinaryString(long num, int width){
		String r=Long.toBinaryString(num);
		while(r.length()<width){r="0"+r;}
		return r;
	}

	/** Formats rows using Arrays.toString, retaining a trailing comma and newline per row. */
	public static final String toString(double[][] a){
		StringBuilder sb=new StringBuilder(256);
		sb.append("[\n");
		for(double[] b : a){
			sb.append(" ").append(Arrays.toString(b)).append(",\n");
		}
		sb.append("]");
		return sb.toString();
	}

	/** Returns the shared reverse map for gene_names_37.1.txt; duplicate names retain one ID. */
	public static final HashMap<String, Integer> geneNameToIdTable(){
		if(geneNameToIdTable==null){
			geneIdToNameTable();
			assert(geneIdToNameTable!=null);
			assert(geneNameToIdTable!=null);
		}
		return geneNameToIdTable;
	}

	/** Lazily reads the fixed legacy gene_names_37.1.txt resource, not the active build's map.
	 * Returns a shared mutable map; complete initialization before publishing to readers. */
	public static final HashMap<Integer, String> geneIdToNameTable(){
		if(geneIdToNameTable==null){

			synchronized(GENEIDLOCK){
				if(geneIdToNameTable==null){

					//			TextFile tf=new TextFile(ROOT_GENE+"gene_names_36.3.txt");
					TextFile tf=new TextFile(ROOT_GENE+"gene_names_37.1.txt", false);
					String[] lines=tf.toStringLines();
					tf.close();
					HashMap<Integer, String> table=new HashMap<Integer, String>((lines.length*3)/2);
					for(String s : lines){
						if(!s.startsWith("#")){
							String[] line=s.split("\t", -1);
							//					assert(line.length==3) : "'"+s+"'";
							if(line.length>=3){

								int key=-1;
								try{
									key=Integer.parseInt(line[1]);
								}catch(NumberFormatException e){
									System.err.println("Bad line: "+s);
									throw new RuntimeException(e);
								}

								table.put(key, (line[2]==null || line[2].length()==0) ? line[1] : line[2]);
							}
						}
					}

					geneIdToNameTable=table;

					HashMap<String, Integer> table2=new HashMap<String, Integer>((lines.length*3)/2);
					for(Integer id : geneIdToNameTable.keySet()){
						table2.put(geneIdToNameTable.get(id), id);
					}
					geneNameToIdTable=table2;
				}
			}
		}
		return geneIdToNameTable;
	}

	/** Lazily returns shared hg18/hg19 chain arrays for build36 to37 or37 to36 only.
	 * Initialization is unsynchronized; other build pairs throw. */
	public static ChainLine[][] getChainLines(int from, int to){
		if(from==36 && to==37){
			if(chains36to37==null){
				chains36to37=ChainBlock.loadChainLines(ROOT_CHAIN+"hg18ToHg19.over.chain");
			}
			return chains36to37;
		}else if(from==37 && to==36){
			if(chains37to36==null){
				chains37to36=ChainBlock.loadChainLines(ROOT_CHAIN+"hg19ToHg18.over.chain");
			}
			return chains37to36;
		}
		throw new RuntimeException("Unknown chain file: "+from+" -> "+to);
	}

	/** Returns "null" or the object's own toString; despite the name, does not traverse arrays. */
	public static final String toStringRecursive(Object a){
		return a==null ? "null" : a.toString();
	}

	/**
	 * Initializes one active genome under the class monitor; the same completed build is a no-op.
	 * Rejects a different build after successful initialization unless state was reset.
	 * A thrown initialization leaves the completion guard unset, allowing another attempt.
	 * Partial public state is not rolled back; ordinary readers remain unsynchronized.
	 * Retries using in-memory metadata must supply fresh complete lists because loading
	 * consumes their entries. Completion means normal loader return; existing read/close
	 * error reporting is unchanged.
	 * @param g Genome build number (must be positive)
	 * @throws RuntimeException if attempting to change already-set genome
	 */
	public static final synchronized void setGenome(int g){
		assert(g>0) : "Genome build number must be at least 1.";
		if(genome_set_to==g){return;}
		if(genome_set_to<0){
			setGenome2(g);
			//DNA-032: Mark completion only after every initialization stage returns normally.
			genome_set_to=g;
		}else{
			throw new RuntimeException("Changing genomes is not currently supported.");
		}
	}

	/**
	 * Internal genome initialization method that loads configuration data.
	 * Reads summary and info files to populate genome statistics and chromosome arrays.
	 * Consumes matching FastaToChromArrays2 metadata lists or reads files from ROOT_GENOME.
	 * GENOME_BUILD is set before I/O for path lookup; the caller marks completion only
	 * after this method returns, including any recursive old-version regeneration.
	 * Exceptions can leave partial public state and consumed metadata lists.
	 * @param g Genome build number to initialize
	 */
	private static final synchronized void setGenome2(int g){
		assert(genome_set_to!=g);
		GENOME_BUILD=g;
		numChroms=-1;
		numBases=-1;
		numDefinedBases=-1;
		numContigs=-1;
		numScaffolds=-1;
		name=null;
		genomeSource=null;
		scaffoldPrefixes=false;
		long fastabytes=-1;
		long fastatime=-1;
		final int currentVersion=FastaToChromArrays2.currentVersion();
		int version=0;

		if(GENOME_BUILD==FastaToChromArrays2.LISTBUILD && FastaToChromArrays2.SUMMARY_LIST!=null){
			for(int i=0; i<FastaToChromArrays2.SUMMARY_LIST.size(); i++){
				final String s=FastaToChromArrays2.SUMMARY_LIST.get(i);
				FastaToChromArrays2.SUMMARY_LIST.set(i, null);
				if(s.charAt(0)=='#'){
					if(s.startsWith("#Version")){
						String[] split=s.split("\t");
						version=(split.length>1 ? Integer.parseInt(split[1]) : 0);
					}
				}else{
					String[] split=s.split("\t");
					String a=split[0];
					String b=split[1];
					if(a.equalsIgnoreCase("chroms")){numChroms=(int)Long.parseLong(b);}
					else if(a.equalsIgnoreCase("bases")){numBases=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("defined")){numDefinedBases=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("contigs")){numContigs=Integer.parseInt(b);}
					else if(a.equalsIgnoreCase("scaffolds")){numScaffolds=Integer.parseInt(b);}
					else if(a.equalsIgnoreCase("interpad")){interScaffoldPadding=Integer.parseInt(b);}
					else if(a.equalsIgnoreCase("undefined")){}
					else if(a.equalsIgnoreCase("name")){name=b;}
					else if(a.equalsIgnoreCase("source")){genomeSource=b;}
					else if(a.equalsIgnoreCase("bytes")){fastabytes=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("last modified")){fastatime=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("scafprefixes")){scaffoldPrefixes=Parse.parseBoolean(b);}
					else{assert(version<currentVersion) : "In array: Unknown term "+s;}
				}
			}
			FastaToChromArrays2.SUMMARY_LIST=null;
		}else{
			String s;
			TextFile tf=new TextFile(ROOT_GENOME+GENOME_BUILD+"/summary.txt", false);
			for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(s.charAt(0)=='#'){
					if(s.startsWith("#Version")){
						String[] split=s.split("\t");
						version=(split.length>1 ? Integer.parseInt(split[1]) : 0);
					}
				}else{
					String[] split=s.split("\t");
					String a=split[0];
					String b=split[1];
					if(a.equalsIgnoreCase("chroms")){numChroms=(int)Long.parseLong(b);}
					else if(a.equalsIgnoreCase("bases")){numBases=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("defined")){numDefinedBases=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("contigs")){numContigs=Integer.parseInt(b);}
					else if(a.equalsIgnoreCase("scaffolds")){numScaffolds=Integer.parseInt(b);}
					else if(a.equalsIgnoreCase("interpad")){interScaffoldPadding=Integer.parseInt(b);}
					else if(a.equalsIgnoreCase("undefined")){}
					else if(a.equalsIgnoreCase("name")){name=b;}
					else if(a.equalsIgnoreCase("source")){genomeSource=b;}
					else if(a.equalsIgnoreCase("bytes")){fastabytes=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("last modified")){fastatime=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("scafprefixes")){scaffoldPrefixes=Parse.parseBoolean(b);}
					else{assert(version<currentVersion) : "In file "+tf.name+": Unknown term "+s;}
				}
			}
			tf.close();
		}
		if(numScaffolds==-1){numScaffolds=numChroms;}

		if(version<currentVersion){
			assert(false) : "The index format has changed in this version of BBTools.  Please delete the /ref/ directory and re-index from the reference fasta," +
					" or use an older version of BBTools.";
			if(new File(ROOT_GENOME+GENOME_BUILD+"/info.txt").exists()){new File(ROOT_GENOME+GENOME_BUILD+"/info.txt").delete();}
			if(new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt").exists()){new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt").delete();}
			if(new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt.gz").exists()){new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt.gz").delete();}
			sysout.println("Regenerating genome info in new format, version "+currentVersion+".");
			dna.FastaToChromArrays2.writeInfo(GENOME_BUILD, numChroms, name, genomeSource, true, scaffoldPrefixes);
			genome_set_to=-1;
			setGenome2(g);
			return;
		}

		assert(numChroms>0 || allowZeroSizedGenome) : "Genome "+g+": numChroms="+numChroms;
		assert(numBases>0 || allowZeroSizedGenome) : "Genome "+g+": numBases="+numBases;
		assert(numDefinedBases>0 || allowZeroSizedGenome) : "Genome "+g+": numDefinedBases="+numDefinedBases;
		assert(numBases>=numDefinedBases) : "Genome "+g+": numBases>numDefinedBases : "+numBases+">"+numDefinedBases;

		chromosomePlusMatrix=new ChromosomeArray[numChroms+1];
		chromLengths=new int[numChroms+1];
		chromDefinedBases=new int[numChroms+1];
		chromUndefinedBases=new int[numChroms+1];
		chromContigs=new int[numChroms+1];
		chromStartPad=new int[numChroms+1];
		chromScaffolds=new int[numChroms+1];

		scaffoldNames=new byte[numChroms+1][][];
		scaffoldLocs=new int[numChroms+1][];
		scaffoldLengths=new int[numChroms+1][];

		if(GENOME_BUILD==FastaToChromArrays2.LISTBUILD && FastaToChromArrays2.INFO_LIST!=null){
			for(int i=0; i<FastaToChromArrays2.INFO_LIST.size(); i++){
				final String s=FastaToChromArrays2.INFO_LIST.get(i);
				FastaToChromArrays2.INFO_LIST.set(i, null);
				if(s.charAt(0)=='#'){
					if(s.startsWith("#Version")){
						String[] split=s.split("\t");
						version=(split.length>1 ? Integer.parseInt(split[1]) : 0);
					}
				}else{
					assert(version==currentVersion);
					String[] split=s.split("\t");
					int chrom=Integer.parseInt(split[0]);
					chromScaffolds[chrom]=Integer.parseInt(split[1]);
					chromContigs[chrom]=(split.length>2 ? Integer.parseInt(split[2]) : -1);
					chromLengths[chrom]=Integer.parseInt(split[3]);
					chromDefinedBases[chrom]=Integer.parseInt(split[4]);
					chromUndefinedBases[chrom]=(split.length>5 ? Integer.parseInt(split[5]) : -1);
					chromStartPad[chrom]=(split.length>6 ? Integer.parseInt(split[6]) : -1);
					//						chromStopPad[chrom]=(split.length>7 ? Integer.parseInt(split[7]) : -1);

				}
			}
			FastaToChromArrays2.INFO_LIST=null;
		}else{
			String s;
			TextFile tf=new TextFile(ROOT_GENOME+GENOME_BUILD+"/info.txt", false);
			for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(s.charAt(0)=='#'){
					if(s.startsWith("#Version")){
						String[] split=s.split("\t");
						version=(split.length>1 ? Integer.parseInt(split[1]) : 0);
					}
				}else{

					if(version>=currentVersion){
						String[] split=s.split("\t");
						int chrom=Integer.parseInt(split[0]);
						chromScaffolds[chrom]=Integer.parseInt(split[1]);
						chromContigs[chrom]=(split.length>2 ? Integer.parseInt(split[2]) : -1);
						chromLengths[chrom]=Integer.parseInt(split[3]);
						chromDefinedBases[chrom]=Integer.parseInt(split[4]);
						chromUndefinedBases[chrom]=(split.length>5 ? Integer.parseInt(split[5]) : -1);
						chromStartPad[chrom]=(split.length>6 ? Integer.parseInt(split[6]) : -1);
//						chromStopPad[chrom]=(split.length>7 ? Integer.parseInt(split[7]) : -1);
					}else{
						tf.close();
						if(new File(ROOT_GENOME+GENOME_BUILD+"/info.txt").exists()){new File(ROOT_GENOME+GENOME_BUILD+"/info.txt").delete();}
						if(new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt").exists()){new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt").delete();}
						if(new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt.gz").exists()){new File(ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt.gz").delete();}
						sysout.println("Regenerating genome info in new format.");
						dna.FastaToChromArrays2.writeInfo(GENOME_BUILD, numChroms, name, genomeSource, true, scaffoldPrefixes);
						tf=new TextFile(ROOT_GENOME+GENOME_BUILD+"/info.txt", false);
					}
				}

			}

			tf.close();
		}

		String fname=ROOT_GENOME+GENOME_BUILD+"/scaffolds.txt.gz";
		boolean hasList=(GENOME_BUILD==FastaToChromArrays2.LISTBUILD && FastaToChromArrays2.SCAF_LIST!=null);

		if(!LOAD_SCAFFOLDS || (!hasList && !new File(fname).exists())){
			for(int i=0; i<scaffoldNames.length; i++){
				scaffoldNames[i]=new byte[][] {("chr"+i).getBytes()};
				scaffoldLocs[i]=new int[] {chromStartPad[i]<0 ? 0 : chromStartPad[i]};
				scaffoldLengths[i]=new int[] {chromLengths[i]};
				if(!LOAD_SCAFFOLDS){chromScaffolds[i]=1;}

				assert(chromScaffolds[i]==1) : "This appears to be an old index version.  " +
					"\nPlease regenerate it from the fasta file by rerunning this program,\nusing the ref=<reference file> and overwrite=true flags.\n"+i+", "+chromScaffolds[i];
			}
		}else{
			for(int chrom=0; chrom<scaffoldNames.length; chrom++){
				int num=chromScaffolds[chrom];
				assert(chrom==0 || num>=1) : chrom+", "+num+", "+Arrays.toString(chromScaffolds);
				if(num>0){
					scaffoldNames[chrom]=new byte[num][];
					scaffoldLocs[chrom]=new int[num];
					scaffoldLengths[chrom]=new int[num];
				}
			}
			int[] count=new int[numChroms+1];

			if(hasList){

				if(verbose){System.err.println("Fetching scaffold names from list:\n\n"+FastaToChromArrays2.SCAF_LIST+"\n\n");}

				for(int i=0; i<FastaToChromArrays2.SCAF_LIST.size(); i++){
					final String s=FastaToChromArrays2.SCAF_LIST.get(i);

					if(verbose){System.err.println("Processing "+s);}

					FastaToChromArrays2.SCAF_LIST.set(i, null);
					if(s.charAt(0)=='#'){
						if(verbose){System.err.println("Control string");}
						if(s.startsWith("#Version")){
							assert(version==currentVersion) : "Wrong index version; please delete /ref/genome/\n"+version+", "+currentVersion;
//							String[] split=s.split("\t");
//							version=(split.length>1 ? Integer.parseInt(split[1]) : 0);
//							assert(version==currentVersion) : "Wrong version: "+version+", "+currentVersion;
						}
					}else{
						String[] split=s.split("\t");
						if(verbose){System.err.println("Split into "+Arrays.toString(split));}
						int chrom=Integer.parseInt(split[0]);
						int x=count[chrom];
						count[chrom]++;

						int scaffoldID=Integer.parseInt(split[1]);
						scaffoldLocs[chrom][x]=Integer.parseInt(split[2]);
						scaffoldLengths[chrom][x]=Integer.parseInt(split[3]);
						scaffoldNames[chrom][x]=split[4].getBytes();
						if(verbose){System.err.println("Set scaffoldNames["+chrom+"]["+x+" to "+(scaffoldNames[chrom][x]==null ? "null" : new String(scaffoldNames[chrom][x])));}
					}
				}
				FastaToChromArrays2.SCAF_LIST=null;
			}else{

				String s;
				TextFile tf=new TextFile(fname, false);
				for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
					if(s.charAt(0)=='#'){
						if(s.startsWith("#Version")){
							assert(version==currentVersion) : "Wrong index version; please delete /ref/genome/\n"+version+", "+currentVersion;
//							String[] split=s.split("\t");
//							version=(split.length>1 ? Integer.parseInt(split[1]) : 0);
//							assert(version==currentVersion) : "Wrong version: "+version+", "+currentVersion;
						}
					}else{
						String[] split=s.split("\t");
						int chrom=Integer.parseInt(split[0]);
						int x=count[chrom];
						count[chrom]++;

						int scaffoldID=Integer.parseInt(split[1]);
						scaffoldLocs[chrom][x]=Integer.parseInt(split[2]);
						scaffoldLengths[chrom][x]=Integer.parseInt(split[3]);
						scaffoldNames[chrom][x]=split[4].getBytes();
					}

				}

				tf.close();
			}
		}

//		assert(false) : (numChroms+1)+", "+(scaffoldLengths==null)+", "+(scaffoldLengths[0]==null)+", "+(scaffoldLengths[1]==null);

//		for(int i=1; i<scaffoldNames.length; i++){
//			stdout.println(Arrays.toString(scaffoldLocs[i]));
//			stdout.println(Arrays.toString(scaffoldLengths[i]));
//			stdout.println(Arrays.toString(scaffoldNames[i]));
//		}

	}

//	public static String contigName(int x){return scaffoldName(x);}
//
//	public static String scaffoldName(int x){
//		if(scaffoldNames==null){return "chr"+x;}
//		return scaffoldNames[x][0];
//	}

	/**
	 * Builds and caches scaffold name lookup table for coordinate translation.
	 * Maps scaffold names to ScafLoc objects containing chromosome and position.
	 * Thread-safe: the fast-path read below relies on scaffoldNameTable being
	 * `volatile` (see its declaration) for a real happens-before guarantee,
	 * AND on this method never publishing the field until it is fully
	 * populated (built in a local variable first) — either alone is
	 * insufficient; both together are the fix (2026-09-02, real bug found
	 * via bbmapnova Phase 0 testing, see the plan for the full writeup).
	 * This publishes the completed map safely, but does not protect later map mutation
	 * or concurrent changes to genome/scaffold metadata. Duplicate names retain the
	 * last encountered location.
	 * @return Shared mutable map, or null when no genome build is active
	 */
	public static HashMap<String, ScafLoc> scafNameTable(){

		if(GENOME_BUILD<0){
			assert(scaffoldNameTable==null);
			return null;
		}
		if(scaffoldNameTable!=null){return scaffoldNameTable;}
		synchronized(SCAFMAPLOCK){
			if(scaffoldNameTable!=null){return scaffoldNameTable;}
			HashMap<String, ScafLoc> local=new HashMap<String, ScafLoc>((int)Tools.min(2L*numScaffolds+10, 1000000000));
			for(int chrom=0; chrom<scaffoldNames.length; chrom++){
				if(scaffoldNames[chrom]!=null){
					for(int scafnum=0; scafnum<scaffoldNames[chrom].length; scafnum++){
						byte[] name=scaffoldNames[chrom][scafnum];
						if(name!=null){
							int loc=scaffoldLocs[chrom][scafnum];
							ScafLoc sc=new ScafLoc(new String(name), chrom, loc);
							local.put(sc.name, sc);
						}
					}
				}
			}
			scaffoldNameTable=local; //Published only now that it is fully populated, into a volatile field.
		}
		return scaffoldNameTable;
	}

	/** Looks up a scaffold using the platform-default byte decoding; null if absent/no build. */
	public static ScafLoc getScafLoc(byte[] name){
		HashMap<String, ScafLoc> map=scafNameTable();
		if(map==null){return null;}
		return map.get(new String(name));
	}
	/** Looks up a scaffold name in the shared map; null if absent or no build is active. */
	public static ScafLoc getScafLoc(String name){
		HashMap<String, ScafLoc> map=scafNameTable();
		if(map==null){return null;}
		return map.get(name);
	}

	/** Returns the borrowed scaffold-name bytes at idx; loc is unused. */
	public static byte[] scaffoldName(int chrom, int loc, int idx){return scaffoldNames[chrom][idx];}
	/** Subtracts the selected scaffold's start; does not validate that loc falls within it. */
	public static int scaffoldRelativeLoc(int chrom, int loc, int idx){return loc-scaffoldLocs[chrom][idx];}

	/**
	 * Finds scaffold index for genomic coordinate within chromosome.
	 * Uses binary search with padding adjustment for gap handling.
	 *
	 * @param chrom Chromosome number
	 * @param loc Genomic coordinate
	 * @return Floor index after adding half the padding, clamped to zero; not proof of containment
	 */
	public static int scaffoldIndex(int chrom, int loc){
		int[] array=scaffoldLocs[chrom];
		if(array==null || array.length<2){return 0;}

		assert(interScaffoldPadding>0);
		loc=loc+interScaffoldPadding/2; //Puts it on closest scaffold if it is between scaffolds

		int idx=Arrays.binarySearch(array, loc);
		if(idx>=0){return idx;} //Perfect hit

		//Otherwise, return closest scaffold.
		int insertPoint=-1-idx;
		assert(insertPoint>=0 && insertPoint<=array.length);
		int r=max(0, insertPoint-1);
		assert(r>=0 && r<array.length);
		assert(r==0 || loc>array[r]);
		assert(r==array.length-1 || loc<array[r+1]);
		return r;
	}

	/** Screens an inclusive interval against the next scaffold boundary with padding tolerance.
	 * Missing/single-scaffold metadata and selection of the final scaffold return true;
	 * this is not a guarantee that every queried base lies inside a scaffold.
	 * TODO: This can be made faster. */
	public static boolean isSingleScaffold(int chrom, int loc1, int loc2){
		assert(loc2>=loc1);
		if(scaffoldLocs==null){return true;}
		assert(chrom>=0 && chrom<scaffoldLocs.length) : chrom+", "+scaffoldLocs.length;
		int[] array=scaffoldLocs[chrom];
		if(array==null || array.length<2){return true;}
		assert(interScaffoldPadding>0);

		int idx=Arrays.binarySearch(array, loc1+interScaffoldPadding);
		final int scaf;
		if(idx>=0){scaf=idx;} //Perfect hit
		else{
			int insertPoint=-1-idx;
			assert(insertPoint>=0 && insertPoint<=array.length);
			scaf=max(0, insertPoint-1);
			assert(scaf>=0 && scaf<array.length);
			assert(scaf==0 || loc1+interScaffoldPadding>array[scaf]);
			assert(scaf==array.length-1 || loc1+interScaffoldPadding<array[scaf+1]);
		}
		if(scaf==array.length-1){return true;}

		int lowerBound=array[scaf]-interScaffoldPadding;
		int upperBound=array[scaf+1];

		if(loc2<lowerBound || loc1>upperBound){return false;} //This could happen if a random read was generated in the start or stop padding.
		assert(scaf==0 || scaf==array.length-1 || (loc1>=lowerBound && loc1<upperBound)) :
			"chrom="+chrom+", loc1="+loc1+", lowerBound="+lowerBound+", loc2="+loc2+", upperBound="+upperBound;
		return loc2<upperBound;
	}

	/** Counts inclusive query bases overlapping one scaffold selected near the query midpoint.
	 * Does not sum overlaps across several scaffolds. Missing/single-scaffold metadata
	 * returns the full query length. Requires loc2>=loc1 and coordinate arithmetic that fits int. */
	public static int scaffoldOverlapLength(int chrom, int loc1, int loc2){
		assert(loc2>=loc1);
		int len=loc2-loc1+1;
		if(scaffoldLocs==null){return len;}
		int[] array=scaffoldLocs[chrom];
		if(array==null || array.length<2){return len;}
		assert(interScaffoldPadding>0);

		int mid=loc1+(interScaffoldPadding+len)/2;
		int idx=Arrays.binarySearch(array, mid);
		final int scaf;
		if(idx>=0){scaf=idx;} //Perfect hit
		else{
			int insertPoint=-1-idx;
			assert(insertPoint>=0 && insertPoint<=array.length);
			scaf=max(0, insertPoint-1);
			assert(scaf>=0 && scaf<array.length);
			assert(scaf==0 || mid>array[scaf]);
			assert(scaf==array.length-1 || mid<array[scaf+1]) : "\nscaf="+scaf+", array.length="+array.length+"\n"+
				"loc1="+loc1+", loc2="+loc2+", mid="+mid+", interScaffoldPadding="+interScaffoldPadding+", "+array[scaf]+", "+array[scaf+1]+"\n"+
				(loc1+interScaffoldPadding)+", "+array[scaf+1];
		}

		int lowerBound=array[scaf];
		int upperBound=lowerBound+scaffoldLengths[chrom][scaf]-1; //overlapLength uses inclusive endpoints
//		assert(upperBound==array[scaf+1]) : lowerBound+", "+upperBound+", "+array[scaf+1]+", "+interScaffoldPadding; //This should fail.

		//A fully contained inclusive query overlaps by its already-computed length.
		if(loc1>=lowerBound && loc2<=upperBound){return len;}
		return Tools.overlapLength(loc1, loc2, lowerBound, upperBound);
	}

	/** Replaces scaffold-name entries with copies truncated at their first whitespace.
	 * Invalidates the name cache only if an entry changes; the next sequential lookup rebuilds it.
	 * Previously returned maps and ScafLoc objects are not updated. Callers must exclude
	 * concurrent lookup and metadata mutation while trimming. */
	public static void trimScaffoldNames(){
		boolean changed=false;
		if(scaffoldNames!=null){
			for(int i=0; i<scaffoldNames.length; i++){
				byte[][] matrix=scaffoldNames[i];
				if(matrix!=null){
					for(int j=0; j<matrix.length; j++){
						byte[] array=matrix[j];
						if(array!=null){
							for(int k=0; k<array.length; k++){
								if(Character.isWhitespace(array[k])){
									matrix[j]=Arrays.copyOf(array, k);
									changed=true;
									break;
								}
							}
						}
					}
				}
			}
		}
		if(changed){scaffoldNameTable=null;}//DNA-033: discard old keys after all names have been replaced.
	}

	/**
	 * Resolves a resource with missing-resource warnings enabled.
	 * @param fname Nonnull filename; only a leading ? requests lookup
	 * @return Resolved path/URL, unchanged literal, or null for an unresolved resource
	 */
	public static final String findPath(String fname){return findPath(fname, true);}

	/**
	 * Locates resource file with optional warning for missing files.
	 * Searches multiple standard locations and JAR resources. Each location tries
	 * the exact name before appending .gz; location priority is unchanged. Literal
	 * paths and names already ending in .gz or .gzip are not expanded.
	 * Network files prefer networks/ before compatible resources/ locations. Local
	 * locations precede classpath/JAR resources, then Dori and NERSC fallbacks.
	 * Literal paths are returned without checking existence; a resolved JAR URL is
	 * not necessarily a filesystem path.
	 *
	 * @param fname Filename to locate (prefix with ? for standard locations)
	 * @param warn Whether to print warnings for missing files
	 * @return Resolved path/URL, unchanged literal, or null for an unresolved resource
	 */
	public static final String findPath(String fname, boolean warn){
		assert(fname!=null);
		if(fname.startsWith("?")){//Look in standard locations
			fname=fname.substring(1);
		}else{//Use this as the literal path
			return fname;
		}
		final boolean network=isNetworkFile(fname);
		final String lower=fname.toLowerCase(Locale.ROOT);
		final boolean tryGzip=!lower.endsWith(".gz") && !lower.endsWith(".gzip");
		String path=(network ? ROOT_NETWORKS : ROOT_RESOURCES)+fname;
		boolean vb=false;
		{
			File f=resourceFile(path, tryGzip);
			path=f.getPath();
			if(!f.exists() && network){
				if(vb){System.err.println("Did not find "+fname+" at "+path);}
				f=resourceFile(ROOT_RESOURCES+fname, tryGzip);path=f.getPath();
			}
			if(!f.exists()){
				if(vb){System.err.println("Did not find "+fname+" at "+path);}
				f=resourceFile(ROOT+fname, tryGzip);
				if(f.exists()){path=f.getPath();}
				else{if(vb){System.err.println("Did not find "+fname+" at "+ROOT);}}
			}
			if(!f.exists()){
				if(vb){System.err.println("Considering fixing "+path+"\n"+path.contains("/file:"));}
				if(path.contains("/file:")){
					String fixed=path.substring(path.lastIndexOf("/file:")+1)+fname;
					f=resourceFile(fixed, tryGzip);
					if(f.exists()){path=f.getPath();}
					else{if(vb){System.err.println("Did not find "+fname+" at "+fixed);}}
				}
			}
			if(!f.exists()){
				if(vb){System.err.println("Considering getResource");}
				URL url=(network ? resourceURL("/networks/"+fname, tryGzip) : null);
				if(url==null){url=resourceURL("/"+fname, tryGzip);}
				if(url!=null){
					String temp=PercentEncoding.codeToSymbol(url.toString());
					if(vb){System.err.println("Found URL "+temp);}
					f=new File(temp);
					path=temp;
				}
			}
			//Central shared resource dirs: local resources/ (above) -> Dori -> NERSC. Each is an absolute
			//path that only exists on its own system, so exists() self-gates it (no host detection). A literal
			//command-line path (no '?') returned at the top and never reaches here, so overrides always win.
			if(!f.exists() && !path.startsWith("jar:")){
				String dori=(network ? DORI_NETWORKS : DORI_RESOURCES)+fname;
				f=resourceFile(dori, tryGzip);
				if(f.exists()){path=f.getPath();}
				else{if(vb){System.err.println("Did not find "+fname+" at "+dori);}}
			}
			if(!f.exists() && network && !path.startsWith("jar:")){
				String dori=DORI_RESOURCES+fname;f=resourceFile(dori, tryGzip);
				if(f.exists()){path=f.getPath();}else{if(vb){System.err.println("Did not find "+fname+" at "+dori);}}
			}
			if(!f.exists() && !path.startsWith("jar:")){
				String nersc=(network ? NERSC_NETWORKS : NERSC_RESOURCES)+fname;
				f=resourceFile(nersc, tryGzip);
				if(f.exists()){path=f.getPath();}
				else{if(vb){System.err.println("Did not find "+fname+" at "+nersc);}}
			}
			if(!f.exists() && network && !path.startsWith("jar:")){
				String nersc=NERSC_RESOURCES+fname;f=resourceFile(nersc, tryGzip);
				if(f.exists()){path=f.getPath();}else{if(vb){System.err.println("Did not find "+fname+" at "+nersc);}}
			}
			if(!f.exists() && !path.startsWith("jar:")){
				if(warn){
					System.err.println("Warning!  Cannot find "+fname+" "+path);
					new Exception().printStackTrace();
				}
				return null;
			}
		}
		if(vb){System.err.println("Found "+fname+" at "+path);}
		return path;
	}

	/** Keeps exact-name preference within one filesystem location. */
	private static File resourceFile(String path, boolean tryGzip){
		assert(path!=null) : "Resource locations must be constructed before testing gzip alternatives";
		final File exact=new File(path);
		if(exact.exists() || !tryGzip){return exact;}
		final File compressed=new File(path+".gz");
		return compressed.exists() ? compressed : exact;
	}

	/** Uses the same exact-then-gzip rule for classpath and JAR resources. */
	private static URL resourceURL(String name, boolean tryGzip){
		assert(name!=null && name.startsWith("/")) : "Classpath resource lookup must be absolute: "+name;
		final URL exact=Primes.class.getResource(name);
		return exact!=null || !tryGzip ? exact : Primes.class.getResource(name+".gz");
	}

	/** True for single- or multi-network files, optionally gzip-compressed. */
	public static boolean isNetworkFile(String fname){
		if(fname==null){return false;}
		final String lower=fname.toLowerCase(Locale.ENGLISH);
		return lower.endsWith(".bbnet") || lower.endsWith(".bbnet.gz") ||
				lower.endsWith(".bbnets") || lower.endsWith(".bbnets.gz");
	}

	/** Returns the smaller int value. */
	public static final int min(int x, int y){return x<y ? x : y;}
	/** Returns the larger int value. */
	public static final int max(int x, int y){return x>y ? x : y;}

	/** Returns the smaller signed byte value. */
	public static final byte min(byte x, byte y){return x<y ? x : y;}
	/** Returns the larger signed byte value. */
	public static final byte max(byte x, byte y){return x>y ? x : y;}

	/** Returns the smaller long value. */
	public static final long min(long x, long y){return x<y ? x : y;}
	/** Returns the larger long value. */
	public static final long max(long x, long y){return x>y ? x : y;}

	/** Returns x if {@code x<y}, otherwise y; unlike Math.min, a NaN in x selects y. */
	public static final double min(double x, double y){return x<y ? x : y;}
	/** Returns x if {@code x>y}, otherwise y; unlike Math.max, a NaN in x selects y. */
	public static final double max(double x, double y){return x>y ? x : y;}

	/** Returns x if {@code x<y}, otherwise y, including comparison ties and NaN. */
	public static final float min(float x, float y){return x<y ? x : y;}
	/** Returns x if {@code x>y}, otherwise y, including comparison ties and NaN. */
	public static final float max(float x, float y){return x>y ? x : y;}

	/** Total number of chromosomes in current genome build */
	public static int numChroms;
	/** Total number of bases across all chromosomes */
	public static long numBases;
	/** Number of defined (non-N) bases across all chromosomes */
	public static long numDefinedBases;
	/** Aggregate contig count from reference metadata. */
	public static int numContigs;
	/** Aggregate scaffold count from reference metadata. */
	public static int numScaffolds;
	/** Nominal N-padding length between packed scaffolds, in bases. */
	public static int interScaffoldPadding;
	/** Per-chromosome length including padding; chromosome number is the array index. */
	public static int[] chromLengths;
	/** Per-chromosome defined-base counts from metadata. */
	public static int[] chromDefinedBases;
	/** Per-chromosome undefined-base counts, including padding. */
	public static int[] chromUndefinedBases;
	/** Per-chromosome contig counts from metadata. */
	public static int[] chromContigs;
	/** Per-chromosome scaffold counts from metadata. */
	public static int[] chromScaffolds;
	/** Per-chromosome leading N-padding lengths. */
	public static int[] chromStartPad;

	/** Allows metadata to describe a genome with zero chromosomes. */
	public static boolean allowZeroSizedGenome=true;

	/** Borrowed name bytes, indexed by chromosome then scaffold. */
	public static byte[][][] scaffoldNames;
	/** Zero-based scaffold starts in packed chromosome coordinates. */
	public static int[][] scaffoldLocs;
	/** Does NOT include interScaffoldPadding */
	public static int[][] scaffoldLengths;
	/** Should be true if scaffold names have extra prefixes (for BBSplitter mode), false otherwise */
	public static boolean scaffoldPrefixes;

	/** Allows translation of sam coordinates back to native coordinates.
	 * volatile is required, not optional (2026-09-02, real bug found via
	 * bbmapnova Phase 0 testing): scafNameTable() below publishes this
	 * field via double-checked locking with an unsynchronized fast-path
	 * read; without volatile, that read has no happens-before guarantee
	 * against the write, and a racing thread can observe a non-null but
	 * not-yet-fully-populated map. Reproduced for real (SamStreamer
	 * workers, "Can't find scaffold in reference with name ..." on an
	 * entry that genuinely exists) before this fix. See
	 * /mnt/c/playground/Nowi/plans/BBMapUpgrade.plan for the full writeup. */
	public static volatile HashMap<String, ScafLoc> scaffoldNameTable;

	/** Reference source text reported by summary metadata; may be null. */
	public static String genomeSource;
	/** Genome name from summary metadata; may be null. */
	public static String name;

	/** Legacy gene-set cache with63 chromosome slots; unloadAll clears every slot.
	 * Range obj1 arrays refer to the shared Gene/GeneSet instances. */
	private static final GeneSet[][] geneSetMatrix=new GeneSet[63][];
	/** Shared legacy genes per chromosome; populated last by loadGenes. */
	private static final Gene[][] geneMatrix=new Gene[63][];
	/** Merged gene-set spans expanded by NEAR, indexed by chromosome. */
	public static final Range[][] geneSetRangeMatrix=new Range[63][];
	/** Merged transcript spans with associated Gene arrays in obj1. */
	public static final Range[][] geneTxRangeMatrix=new Range[63][];
	/** Legacy merged coding spans; ordering can be invalid as recorded in DNA-036. */
	public static final Range[][] geneCodeRangeMatrix=new Range[63][];
	/** Merged coding/exon intersections with associated Gene arrays in obj1. */
	private static final Range[][] geneCodeAndExonRangeMatrix=new Range[63][];
	/** Merged whole-exon spans with associated Gene arrays in obj1. */
	public static final Range[][] exonRangeMatrix=new Range[63][];
	/** Shared active-build ID-to-gene-set map, invalidated by unloadAll.
	 * Nonvolatile despite the double-checked loader: initialize during setup before
	 * publishing to workers, as required for chromosomePlusMatrix below. */
	public static HashMap<Integer, ArrayList<GeneSet>> geneIDTable;

	/** Ranges within genes and exons or within NEAR their ends */
	public static final Range[][] geneNearbyRangeMatrix=new Range[63][];

	/** Lazily-loaded chromosome arrays, indexed by chromosome number. The double-checked-locking
	 * accessors (getChromosome) read this on the fast path without volatile/sync; safe because
	 * chromosomes are loaded during single-threaded setup before worker threads fork (happens-before).
	 * Not safe to FIRST-load a chromosome concurrently from multiple already-running threads. */
	public static ChromosomeArray[] chromosomePlusMatrix;

	/** Shared ID-to-name map from gene_names_37.1.txt, independent of active build. */
	private static HashMap<Integer, String> geneIdToNameTable;
	/** Reverse legacy name map; duplicate names retain one ID. */
	private static HashMap<String, Integer> geneNameToIdTable;

	/** Per-chromosome symbol/accession maps; unloadGenes does not clear them. */
	private static final HashMap<String, GeneSet>[] geneNameTable=new HashMap[63];
	/** Per-chromosome transcript-accession maps; missing accessions use chromosome_id. */
	private static final HashMap<String, Gene>[] transcriptNameTable=new HashMap[63];

	/** Lazily loaded hg18-to-hg19 chain data, shared with callers. */
	public static ChainLine[][] chains36to37;
	/** Lazily loaded hg19-to-hg18 chain data, shared with callers. */
	public static ChainLine[][] chains37to36;

	/** Selects inclusive transcript coordinates in findGeneRanges. */
	private static final int TX_RANGE=0;
	/** Selects inclusive coding coordinates in findGeneRanges. */
	private static final int CODE_RANGE=1;

	/** Distance in bases for "nearby" gene region calculations */
	public static final int NEAR=200;

	/** Returns the installation root determined once from the Data class location. */
	public static final String ROOT(){return ROOT;}
	/** Returns the installation's resources directory prefix. */
	public static final String RESOURCES(){return ROOT_RESOURCES;}
	/** Returns the installation's networks directory prefix. */
	public static final String NETWORKS(){return ROOT_NETWORKS;}

	/** Found dynamically.  Path to /BBTools/ (the root of the project, NOT including /current/) */
	private static final String ROOT;
	/** Path to /BBTools/resources/ */
	private static final String ROOT_RESOURCES;
	/** Path to /BBTools/networks/ */
	private static final String ROOT_NETWORKS;

	/** Central shared resource directory on JGI's Dori cluster, checked after local resources/ (findPath).
	 *  Hardcoded default; an absolute path that only exists on Dori, so it self-gates via exists(). */
	public static String DORI_RESOURCES="/clusterfs/jgi/groups/gentech/homes/bbushnell/resources/";
	/** Dori network-file fallback prefix, used only when the candidate exists. */
	public static String DORI_NETWORKS="/clusterfs/jgi/groups/gentech/homes/bbushnell/networks/";
	/** Central shared resource directory on NERSC, checked last (after Dori). Hardcoded default;
	 *  self-gates via exists() (the path only exists on NERSC). */
	public static String NERSC_RESOURCES="/global/cfs/cdirs/bbtools/resources/";
	/** NERSC network-file fallback prefix, checked after local/classpath and Dori locations. */
	public static String NERSC_NETWORKS="/global/cfs/cdirs/bbtools/networks/";

	/** Mutable reference base configured by setPath; separate
	 * from the immutable installation/resource roots above. */
	public static String ROOT_BASE;
	/** Reference-data directory prefix, normally ROOT_BASE plus ref/. */
	public static String ROOT_REF;
	/** Genome-build directory prefix beneath the reference-data root. */
	public static String ROOT_GENOME;
	/** Index directory prefix; literal null mode leaves its initial ref/index/ value. */
	public static String ROOT_INDEX;
	/** Legacy gene-annotation directory prefix. */
	public static String ROOT_GENE;
	/** Legacy genome-chain directory prefix. */
	public static String ROOT_CHAIN;
	/** Temporary-directory prefix; custom reference prefixes retain the previous value. */
	public static String ROOT_TEMPDIR;
	/** Working directory captured from user.dir when setPath runs. */
	public static String ROOT_CURRENT;
	/** Reference quality-data directory prefix. */
	public static String ROOT_QUALITY;

	static{
		String s=new File(Data.class.getClassLoader().getResource(
			Data.class.getName().replace('.', '/') + ".class").getFile()).getAbsolutePath();
		s=PercentEncoding.codeToSymbol(s);
		s=s.replace('\\', '/');
		// Fix garbled path when running from jar: getFile() returns "file:/path/jar!/class"
		// new File() prepends cwd, giving "cwd/file:/path/jar!/class"
		if(s.contains("/file:")){
			s=s.substring(s.lastIndexOf("/file:")+6); // strips to "/path/bbtools.jar!/..."
		}
		String root=s.replace("dna/Data.class", "");
		if(root.contains("bbtools.jar!")){root=root.substring(0, root.indexOf("bbtools.jar!"));}
		if(root.endsWith("/current/") || root.endsWith("/current")){
			root=root.substring(0, root.lastIndexOf("/current"));
		}
		if(!root.endsWith("/")){root=root+"/";}
		// ROOT now = "/path/bbtools/bbtools.jar!/" → parent = "/path/bbtools/" → resources found! ✓
		setPath(Shared.WINDOWS ? "?Shared.WINDOWS" : "?unix");
		if(!Shared.WINDOWS || true){setPath("?local");}
		ROOT=root;
		ROOT_RESOURCES=ROOT+"resources/";
		ROOT_NETWORKS=ROOT+"networks/";
	}

	/**
	 * Configures reference paths; immutable installation/resource roots are unchanged.
	 * Local aliases use relative ref/ paths; ?unix additionally consults TEMPDIR.
	 * ?Shared.WINDOWS keeps relative reference roots and selects C:/workspace/tempdir/.
	 * Custom prefixes replace reference roots but retain the previous temporary path.
	 * Null or literal "null" clears most roots but leaves ROOT_INDEX at ref/index/.
	 * An empty string remains a custom / prefix. Selector comparisons use the default locale.
	 * @param path Null, ?local, ?unix, ?Shared.WINDOWS, local alias, literal "null", or custom prefix
	 */
	public static void setPath(String path){
//		System.err.println("***"+path);
		//DNA-026: Preserve typed null until the null-mode selection below.
		if(path!=null && path.indexOf('\\')>=0){path=path.replace('\\', '/');}
		String mode=(path==null ? "null" : path.toLowerCase());
		boolean local=mode.equals("?local") || mode.equals(".") || mode.equals("/.") || mode.equals("./");
		//DNA-030: Compare the already-lowercased mode to the lowercase selector.
		boolean win=mode.equals("?shared.windows");
		boolean unix=mode.equals("?unix");

		ROOT_CURRENT=System.getProperty("user.dir");

		ROOT_BASE="";
		ROOT_REF="ref/";
		ROOT_GENOME=ROOT_REF+"genome/";
		ROOT_INDEX=ROOT_REF+"index/";
		ROOT_GENE=ROOT_REF+"genes/";
		ROOT_CHAIN=ROOT_REF+"chain/";
		ROOT_QUALITY=ROOT_REF+"qual/";

		if(local){
			ROOT_TEMPDIR=ROOT_BASE;
		}else if(win){
			ROOT_TEMPDIR="C:/workspace/tempdir/";
		}else if(unix){
			String s=System.getenv().get("TEMPDIR");
			ROOT_TEMPDIR=(s==null ? ROOT_BASE : s+"/");
		}else if(!"null".equals(mode)){
			if(!path.endsWith("/")){path=path+"/";}
			ROOT_BASE=path;
			ROOT_REF=path+"ref/";
			ROOT_GENOME=ROOT_REF+"genome/";
			ROOT_INDEX=ROOT_REF+"index/";
			ROOT_GENE=ROOT_REF+"genes/";
			ROOT_CHAIN=ROOT_REF+"chain/";
			ROOT_QUALITY=ROOT_REF+"qual/";
		}else{
			ROOT_BASE=null;
			ROOT_REF=null;
			ROOT_GENOME=null;
			ROOT_GENE=null;
			ROOT_CHAIN=null;
			ROOT_QUALITY=null;
			ROOT_TEMPDIR=null;
		}
	}

	/** Currently active genome build number (-1 if unset) */
	public static int GENOME_BUILD=-1;
	/** Completed build guard, assigned after normal initialization return (DNA-032); -1 permits an attempt. */
	private static int genome_set_to=-1;

	/** Compile-time diagnostic switch. */
	public static final boolean verbose=false;

	/** seqGene, knownGene, refGene, unionGene, seqRefGene, ccs */
	public static String GENE_MAP="seqRefGene";

	/** Serializes construction of the legacy gene-ID and name maps. */
	private static final String GENEIDLOCK=new String("GENEIDLOCK");

	/** Striped per-chromosome loader monitors; the array itself is the separate bulk-load monitor. */
	private static final String[] CHROMLOCKS=new String[256];

	static{
		for(int i=0; i<CHROMLOCKS.length; i++){
			CHROMLOCKS[i]=new String(i+"");
		}
	}

	/** Initial string-cache capacity; literal avoids circular initialization through Primes. */
	private static final int INTERN_MAP_SIZE=20011; //Do NOT look up in primes; causes an init ordering issue.
	/** Soft entry-count threshold triggering a cache clear; concurrent inserts can overshoot. */
	private static final int INTERN_MAP_LIMIT=(int)(INTERN_MAP_SIZE*0.75f);

	/** Shared equality cache for short strings; clearing ends its reference-identity guarantee. */
	private static final ConcurrentHashMap<String, String> INTERNMAP=new ConcurrentHashMap<String, String>(INTERN_MAP_SIZE);
//	public static final void unloadInternMap(){
//		INTERNMAP=new HashMap<String, String>(INTERN_MAP_SIZE);
//	}

	/** Interns all strings in array to reduce memory usage.
	 * @param s Array of strings to intern (modified in-place) */
	public static final void intern(String[] s){
		if(s==null){return;}
		for(int i=0; i<s.length; i++){s[i]=intern(s[i]);}
	}
	/**
	 * Shares selected short strings through a bounded, periodically cleared cache.
	 * Null and lengths above25 are unchanged. Base-only strings of length7..25
	 * receive a new String; other eligible strings use forceIntern.
	 *
	 * @param s String to intern
	 * @return Equal-content string; reference identity is not guaranteed across cache clears
	 */
	public static String intern(String s){
		if(s==null || s.length()>25){return s;}
//		calls++;
//
//		if(s.length()>0 && s.charAt(0)!='?'){
//			s=condense(s);
//		}

		if(s.length()<2){
//			return s.intern();
			return forceIntern(s);
		}
		boolean acgtn=AminoAcid.containsOnlyACGTNQ(s);

		if(acgtn){
			if(s.length()<4){
//				return s.intern();
				return forceIntern(s);
			}
			if(s.length()>6){
				return new String(s);
			}
		}

		//Otherwise it is a non-base string of length2..25, or a base string of length4..6.
		return forceIntern(s);
	}

	/**
	 * Shares a nonnull string through the concurrent map, clearing it above a soft size limit.
	 * Concurrent insertion can overshoot the limit, and clearing ends any cache-wide
	 * identity guarantee. The diagnostic calls counter is not atomic.
	 *
	 * @param s0 String to intern
	 * @return Existing equal string in the current cache, or s0 when newly inserted
	 */
	public static String forceIntern(final String s0){
		assert(s0!=null);
		calls++;

//		if(s.length()<2){return s.intern();}
//		boolean acgtn=AminoAcid.containsOnlyACGTNQ(s);
//
//		if(acgtn){
//			if(s.length()<4){return s.intern();}
//		}

		if(INTERNMAP.size()>INTERN_MAP_LIMIT){
			synchronized(INTERNMAP){
				if(INTERNMAP.size()>INTERN_MAP_LIMIT){
//					System.err.println("INTERNMAP overflow caused by "+s0);
					INTERNMAP.clear();
				}
			}
		}
		String s=INTERNMAP.putIfAbsent(s0, s0);
		return s==null ? s0 : s;
	}
	/** Nonatomic diagnostic count of forceIntern calls; not a reliable concurrent total. */
	static int calls=0;

	/** System output stream for status messages */
	public static PrintStream sysout=System.err;//System.out;

	/** Whether chromosome files are gzip compressed */
	public static boolean CHROMGZ=true;
	/** Whether to load detailed scaffold information from files */
	public static boolean LOAD_SCAFFOLDS=true;

//	private static final boolean GUNZIP=testExecute("gunzip --help");
//	private static final boolean GZIP=testExecute("gzip --help");
//	private static final boolean SAMTOOLS=testExecute("samtools --help");

	/** Returns the cached gunzip setting; disabled initially, with zero delegating to GZIP. */
	public static boolean GUNZIP(){return GUNZIP==0 ? GZIP() : GUNZIP>0;}
//	public static boolean UNPIGZ(){return UNPIGZ==0 ? PIGZ() : UNPIGZ>0;}

	/** Probes bgzip once outside Windows and parses version1.x option heuristics.
	 * Availability does not guarantee every version-specific option is recognized.
	 * Note:
	 * Magic number of bgzip files is (first 4 bytes):
	 * 1f 8b 08 04
	 * 31 139 8 4
	 *  = 529205252
	 *
	 * gzip/pigz:
	 * 1f 8b 08 00 (actually I've seen it as 1f 8b 08 08 also)
	 * 31 139 8 0
	 *  = 529205248
	 *
	 * od --format=x1 --read-bytes=16 names.txt_gzip.gz
	 */
	public static boolean BGZIP(){
//		System.err.println("Looking for bgzip.");
		if(BGZIP==0 && !Shared.WINDOWS){
			synchronized(SUBPROCSYNC){
				if(BGZIP==0){
					ByteBuilder bb=new ByteBuilder();
					BGZIP=testExecute("bgzip -h", bb);

					try{
						if(BGZIP>0){
							String[] lines=bb.toString().split("\n");
							for(String line : lines){
								if(line.startsWith("Version:")){
									String[] split=line.split("\\s+");
									if(split.length>1){
										BGZIP_VERSION=split[1];

										BGZIP_VERSION_threadsFlag=false;
										BGZIP_VERSION_levelFlag=false;

										String[] ss=BGZIP_VERSION.split("[-\\.]");

//										System.err.println(Arrays.toString("1.7-8".split("[-\\.]")));

										int[] is=new int[3];
										for(int i=0; i<3 && i<ss.length; i++){
											String s=ss[i];
											int x=Integer.parseInt(s);
											is[i]=x;
										}

										if(is[0]<1){
											//do nothing
										}else if(is[1]>=9){
											BGZIP_VERSION_threadsFlag=BGZIP_VERSION_levelFlag=true;
										}else if(is[1]>=4){
											BGZIP_VERSION_threadsFlag=true;
										}else if(is[1]>=3 && is[2]>=1){
											BGZIP_VERSION_threadsFlag=true;
										}
									}
								}
							}
//							System.err.println("Found BGZIP"+(BGZIP_VERSION==null ? "." : " "+BGZIP_VERSION));
						}else{
//							System.err.println("Could not find BGZIP.");
						}
					}catch(NumberFormatException e){
						System.err.println("Warning - trouble parsing BGZIP version:\n"+BGZIP_VERSION);
					}

				}
			}
		}
//		System.err.println("BGZIP="+BGZIP);
//		assert(BGZIP>0) : "bgzip not found.";
		return BGZIP>0;
	}
	/** Requires bgzip availability and the parsed threads-option capability. */
	public static boolean BGZIP_THREADED(){
		return BGZIP() && BGZIP_VERSION_threadsFlag;
	}
	/** Returns cached gzip availability, probing --version outside Windows if unknown. */
	public static boolean GZIP(){
		if(GZIP==0 && !Shared.WINDOWS){
			synchronized(SUBPROCSYNC){
				if(GZIP==0){
					ByteBuilder bb=new ByteBuilder();
					GZIP=testExecute("gzip --version", bb);
				}
			}
		}
		return GZIP>0;
	}
	/** Probes pigz --version once, including on Windows, and records parsed option capabilities. */
	public static boolean PIGZ(){
		if(PIGZ==0){
			synchronized(SUBPROCSYNC){
				if(PIGZ==0){
					ByteBuilder bb=new ByteBuilder();
					PIGZ=testExecute("pigz --version", bb);

					try{
						if(PIGZ>0){
							String[] lines=bb.toString().split("\n");
							for(String line : lines){
								if(line.startsWith("pigz")){
									String[] split=line.split("\\s+");
									if(split.length>1){
										PIGZ_VERSION=split[1];
										PIGZ_VERSION_23plus=false;
										PIGZ_VERSION_231plus=false;

										String[] ss=PIGZ_VERSION.split("\\.");
										int[] is=new int[3];
										for(int i=0; i<3 && i<ss.length; i++){
											String s=ss[i];
											int x=Integer.parseInt(s);
											is[i]=x;
										}

										if(is[0]>2){
											PIGZ_VERSION_231plus=PIGZ_VERSION_23plus=true;
										}else if(is[0]==2){
											if(is[1]>3){
												PIGZ_VERSION_231plus=PIGZ_VERSION_23plus=true;
											}else if(is[1]==3){
												PIGZ_VERSION_23plus=true;
												if(is[2]>=1){
													PIGZ_VERSION_231plus=true;
												}
											}
										}
									}
								}
							}
//							System.err.println("Found pigz"+(PIGZ_VERSION==null ? "." : " "+PIGZ_VERSION));
						}else{
//							System.err.println("Could not find pigz.");
						}
					}catch(NumberFormatException e){
						System.err.println("Warning - trouble parsing pigz version:\n"+PIGZ_VERSION);
					}
				}
			}
		}
		return PIGZ>0;
	}
	/** Returns cached dsrc availability, probing --version if unknown, including on Windows. */
	public static boolean DSRC(){
		if(DSRC==0){
			synchronized(SUBPROCSYNC){
				if(DSRC==0){
					ByteBuilder bb=new ByteBuilder();
					DSRC=testExecute("dsrc --version", bb);
				}
			}
		}
		return DSRC>0;
	}
	/** Returns cached bzip2 availability, probing --version outside Windows if unknown. */
	public static boolean BZIP2(){
		if(BZIP2==0 && !Shared.WINDOWS){
			synchronized(SUBPROCSYNC){
				if(BZIP2==0){
					ByteBuilder bb=new ByteBuilder();
					BZIP2=testExecute("bzip2 --version", bb);
				}
			}
		}
		return BZIP2>0;
	}
	/** Returns cached pbzip2 availability, probing --version outside Windows if unknown. */
	public static boolean PBZIP2(){
		if(PBZIP2==0 && !Shared.WINDOWS){
			synchronized(SUBPROCSYNC){
				if(PBZIP2==0){
					ByteBuilder bb=new ByteBuilder();
					PBZIP2=testExecute("pbzip2 --version", bb);
				}
			}
		}
		return PBZIP2>0;
	}
	/** Returns cached lbzip2 availability, probing --version outside Windows if unknown. */
	public static boolean LBZIP2(){
		if(LBZIP2==0 && !Shared.WINDOWS){
			synchronized(SUBPROCSYNC){
				if(LBZIP2==0){
					ByteBuilder bb=new ByteBuilder();
					LBZIP2=testExecute("lbzip2 --version", bb);
				}
			}
		}
		return LBZIP2>0;
	}
	/** Whether native BAM output is allowed or samtools is available. */
	public static boolean BAM_SUPPORT_OUT(){
		return ReadWrite.ALLOW_NATIVE_BAM_OUT || SAMTOOLS();
	}
	/** Whether native BAM input is allowed or samtools/sambamba is available. */
	public static boolean BAM_SUPPORT_IN(){
		return ReadWrite.ALLOW_NATIVE_BAM_IN || SAMTOOLS() || SAMBAMBA();
	}
	/** Probes samtools once outside Windows, honoring USE_SAMTOOLS before the first probe. */
	public static boolean SAMTOOLS(){
		if(SAMTOOLS==0 && !Shared.WINDOWS){
			synchronized(SUBPROCSYNC){
				if(!USE_SAMTOOLS){SAMTOOLS=-1;}
				if(SAMTOOLS==0){
					ByteBuilder bb=new ByteBuilder();
					SAMTOOLS=testExecute("samtools", bb);

					if(SAMTOOLS>0){
						String[] lines=bb.toString().split("\n");
						for(String line : lines){
							if(line.startsWith("Version:")){
								String[] split=line.split("\\s+");
								if(split.length>1){
									SAMTOOLS_VERSION=split[1];
									SAMTOOLS_VERSION_1x=SAMTOOLS_VERSION.startsWith("1.");
								}
							}
						}
						System.err.println("Found samtools"+(SAMTOOLS_VERSION==null ? "." : " "+SAMTOOLS_VERSION));
					}else{
						System.err.println("Could not find samtools.");
					}
				}
			}
		}
		return SAMTOOLS>0;
	}
	/** Probes sambamba once outside Windows, honoring USE_SAMBAMBA before the first probe. */
	public static boolean SAMBAMBA(){
		if(SAMBAMBA==0 && !Shared.WINDOWS){
			synchronized(SUBPROCSYNC){
				if(!USE_SAMBAMBA){SAMBAMBA=-1;}
				if(SAMBAMBA==0){
					ByteBuilder bb=new ByteBuilder();
					SAMBAMBA=testExecute("sambamba -q", bb);
					if(SAMBAMBA>0){System.err.println("Found sambamba.");}else{System.err.println("Could not find sambamba.");}
				}
			}
		}
		return SAMBAMBA>0;
	}
	/** Returns cached sh availability, probing --version outside Windows if unknown. */
	public static boolean SH(){
		if(SH==0 && !Shared.WINDOWS){
			synchronized(SUBPROCSYNC){
				if(SH==0){
					ByteBuilder bb=new ByteBuilder();
					SH=testExecute("sh --version", bb);
				}
			}
//			System.err.println(SH>0 ? "Found sh." : "Could not find sh.");
			if(SH<0){System.err.println("Could not find sh; won't launch I/O subprocesses.");}
		}
		return SH>0;
	}
	/** Serializes external-tool availability probes and associated version parsing. */
	private static final String SUBPROCSYNC=new String("SUBPROCSYNC");
	/** Serializes scaffold-name map construction before volatile publication. */
	private static final String SCAFMAPLOCK=new String("SCAFMAPLOCK");

	/** Gunzip setting, disabled initially. Probe states: 0=not checked, -1=unavailable/disabled, 1=available.
	 * Initialize before sharing with workers; fields are not volatile and public flags
	 * do not invalidate completed probes. Set these to zero to enable a new probe. */
	private static int GUNZIP=-1;
//	private static int UNPIGZ=0;
	/** Cached gzip availability: 0 unknown, -1 unavailable, 1 available. */
	private static int GZIP=0;
	/** Cached bgzip availability: 0 unknown, -1 unavailable, 1 available. */
	private static int BGZIP=0;
	/** Cached pigz availability: 0 unknown, -1 unavailable, 1 available. */
	private static int PIGZ=0;
	/** Cached dsrc availability: 0 unknown, -1 unavailable, 1 available. */
	private static int DSRC=0;
	/** Cached bzip2 availability: 0 unknown, -1 unavailable, 1 available. */
	private static int BZIP2=0;
	/** Cached pbzip2 availability: 0 unknown, -1 unavailable, 1 available. */
	private static int PBZIP2=0;
	/** Cached lbzip2 availability: 0 unknown, -1 unavailable, 1 available. */
	private static int LBZIP2=0;
	/** Cached samtools availability: 0 unknown, -1 unavailable/disabled, 1 available. */
	private static int SAMTOOLS=0;
	/** Cached sambamba availability: 0 unknown, -1 unavailable/disabled, 1 available. */
	private static int SAMBAMBA=0;
	/** Allows samtools at its first probe; changing this does not invalidate a cached result. */
	public static boolean USE_SAMTOOLS=true;
	/** Allows sambamba at its first probe; changing this does not invalidate a cached result. */
	public static boolean USE_SAMBAMBA=true;
	/** Parsed samtools version text, or null before a successful version parse. */
	public static String SAMTOOLS_VERSION=null;
	/** Legacy samtools format selector; defaults true until a parsed version says otherwise. */
	public static boolean SAMTOOLS_VERSION_1x=true;
	/** Parsed pigz version text, or null before a successful version parse. */
	public static String PIGZ_VERSION=null;
	/** Parsed pigz capability for the iteration-count option. */
	public static boolean PIGZ_VERSION_231plus=false; //Allows -I flag
	/** Parsed pigz capability for compression level11. */
	public static boolean PIGZ_VERSION_23plus=false; //Allows -11 flag
	/** Parsed bgzip version text, or null before a successful version parse. */
	public static String BGZIP_VERSION=null;
	/** Parsed bgzip capability for the thread-count option; see BGZIP's version heuristic. */
	public static boolean BGZIP_VERSION_threadsFlag=false; //Allows -@ flag; not sure about the minimum version, but 1.4 is capable
	/** Parsed bgzip capability for the compression-level option. */
	public static boolean BGZIP_VERSION_levelFlag=false; //Allows -l flag; not sure about the minimum version, but 1.9 is capable
	/** Cached shell availability: 0 unknown, -1 unavailable, 1 available. */
	private static int SH=0;

	/** Unused legacy probe: consumes stderr only, so a stdout-heavy child can block. */
	private static int testExecute_old(String s, ByteBuilder bb){
//		System.err.println("Testing "+s);
		try{
			Process p=Runtime.getRuntime().exec(s);
//			System.err.println("Got process.");

			InputStream stream=p.getErrorStream();
			if(bb==null){
				while(stream.read()>-1){}
			}else{
				bb.clear();
				for(int b=stream.read(); b>=0; b=stream.read()){
					bb.append((char)b);
//					System.err.print((char)b);
//					System.err.print(""+(int)b);
				}
			}

			//TODO: Probable bug - DNA-035 (historical DNA-029): an already-exited child skips waitFor and its real exit status is lost.
			int exitValue=0;
			while(p.isAlive()){
				try{
					exitValue=p.waitFor();
				}catch(InterruptedException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			if(exitValue>=126){return -1;}

//			return p.exitValue()==0;
//			System.err.println("This system does has "+s+" installed.");
		}catch(IOException e){
//			System.err.println("This system does not have "+s+" installed.");
			// TODO Auto-generated catch block
//			e.printStackTrace();
			return -1;
		}
		return 1;
	}

	/** Probes a whitespace-separated command without shell quoting/expansion.
	 * Georg Rath version: closes stdin and consumes merged stdout/stderr to avoid
	 * child pipe blockage. If supplied, bb is cleared and receives the merged text.
	 * Returns -1 on launch/read failure or exit>=126, otherwise1;
	 * help commands may return nonzero. Always collects the terminal exit status. */
	private static int testExecute(String s, ByteBuilder bb){
//		System.err.println("Testing "+s);
		try{
			ProcessBuilder pb=new ProcessBuilder(s.split("\\s+"));
//			System.err.println("Made builder.");
			pb.redirectErrorStream(true); //redirect to stdout
//			System.err.println("Redirected stderr.");
			Process p=pb.start();
//			System.err.println("Got process.");
			p.getOutputStream().close(); //close stdin
//			System.err.println("Closed stdin.");
			//			InputStream stream=p.getErrorStream();
			InputStream stream=p.getInputStream(); //stdout
//			System.err.println("Grabbed stdout.");

			if(bb==null){
//				System.err.println("Reading to null.");
				while(stream.read()>-1){}
			}else{
//				System.err.println("Reading to buffer.");
				bb.clear();
				for(int b=stream.read(); b>=0; b=stream.read()){
					bb.append((char)b);
					// System.err.print((char)b);
					// System.err.print(""+(int)b);
				}
			}
//			System.err.println("Finished reading.");
			int exitValue=0;
			while(p.isAlive()){
				try{
					exitValue=p.waitFor();
				}catch(InterruptedException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
//			System.err.println("Process terminated.");
			//The child may have exited while output drained, bypassing the wait loop (DNA-029).
			exitValue=p.exitValue();
			if(exitValue>=126){
//				System.err.println("Returning -1 due to 126.");
				return -1;
			}
			// return p.exitValue()==0;
			// System.err.println("This system has "+s+" installed.");
		}catch(IOException e){
//			System.err.println("Returning -1 due to exception: "+e);
			// System.err.println("This system does not have "+s+" installed.");
			// TODO Auto-generated catch block
			// e.printStackTrace();
			return -1;
		}
//		System.err.println("Returning 1");
		return 1;
	}

}
