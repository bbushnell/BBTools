package ddl;

import java.io.PrintStream;
import java.util.ArrayList;

import cardinality.DynamicDemiLog;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.IntHashSet;
import map.LongHashSet;
import map.LongObjectMap;
import parse.Parse;
import parse.Parser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import tax.TaxNode;
import tax.TaxTree;

/**
 * Builds a DDL kmer blacklist from pre-built DDL sketch files.
 * Reads sketches with kmer arrays (kmers=t), promotes each TaxID
 * to genus level, counts distinct genera per kmer, and outputs
 * kmers exceeding the threshold as FASTA with multi-level headers.
 *
 * Accepts multiple comma-separated input files; processes each
 * sequentially and discards sketch data after kmer extraction
 * to minimize memory usage.
 *
 * Condense mode (condense=t blacklist=X): condenses double-sized
 * sketches to half size while preferring non-blacklisted kmers,
 * then validates. Approximates full re-sketching without reading
 * raw genome data.
 *
 * @author Brian Bushnell, Noire
 * @date May 24, 2026
 */
public class DDLBlacklistMaker {

	public static void main(String[] args){
		Timer t=new Timer();
		DDLBlacklistMaker bm=new DDLBlacklistMaker(args);
		bm.process(t);
	}

	public DDLBlacklistMaker(String[] args){
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("mintaxcount") || a.equals("mintaxa") || a.equals("mincount")){
				minTaxCount=Integer.parseInt(b);
			}else if(a.equals("minentropy") || a.equals("entropy")){
				minEntropy=Float.parseFloat(b);
			}else if(a.equals("entropyk") || a.equals("ek")){
				entropyK=Integer.parseInt(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("validate")){
				validate=Parse.parseBoolean(b);
			}else if(a.equals("samples") || a.equals("validatesamples")){
				validateSamples=Integer.parseInt(b);
			}else if(a.equals("exponent") || a.equals("ebits")){
				DynamicDemiLog.setExponent(Integer.parseInt(b));
			}else if(a.equals("blacklist") || a.equals("bl")){
				blacklistFile=b;
			}else if(a.equals("condense")){
				condense=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		{
			Parser.processQuality();
			overwrite=parser.overwrite;
			out=parser.out1;
			in=parser.in1;
		}

		assert(in!=null) : "No input file specified.";
		assert(out!=null || validate || condense) : "No output file or action specified.";
	}

	void process(Timer t){
		if(condense){
			processCondense(t);
		}else{
			processBlacklist(t);
		}
	}

	/*--------------------------------------------------------------*/
	/*-----------        Condense + Validate        ----------------*/
	/*--------------------------------------------------------------*/

	void processCondense(Timer t){
		final LongHashSet bl;
		if(blacklistFile!=null){
			DynamicDemiLog.loadBlacklist(blacklistFile);
			bl=getBlacklistSet();
		}else{
			bl=null;
		}

		final String[] inFiles=in.split(",");
		final ArrayList<DDLRecord> condensed=new ArrayList<DDLRecord>();
		long totalRecords=0;
		int blacklisted=0, kept=0;

		for(String inFile : inFiles){
			outstream.println("Loading "+inFile.trim());
			ArrayList<DDLRecord> records=DDLLoader.loadFile(inFile.trim(), k);
			outstream.println("  "+records.size()+" records");
			totalRecords+=records.size();

			for(DDLRecord rec : records){
				DynamicDemiLog ddl=rec.ddl;
				int oldBuckets=ddl.maxArray().length;
				if(oldBuckets<2){continue;}
				int newBuckets=oldBuckets/2;

				char[] oldMax=ddl.maxArray();
				long[] oldKmers=ddl.kmerArray();

				char[] newMax=new char[newBuckets];
				for(int i=0; i<newBuckets; i++){
					int j=i+newBuckets;
					char scoreA=oldMax[i], scoreB=oldMax[j];
					boolean blA=(bl!=null && oldKmers!=null && oldKmers[i]!=0 && bl.contains(oldKmers[i]));
					boolean blB=(bl!=null && oldKmers!=null && oldKmers[j]!=0 && bl.contains(oldKmers[j]));

					if(blA && !blB){
						newMax[i]=scoreB; blacklisted++;
					}else if(blB && !blA){
						newMax[i]=scoreA; blacklisted++;
					}else{
						newMax[i]=(scoreA>=scoreB ? scoreA : scoreB);
						if(blA && blB){kept++;}
					}
				}

				DynamicDemiLog cddl=DynamicDemiLog.fromArray(newMax, (int)rec.id, rec.name, k, -1);
				condensed.add(new DDLRecord(cddl, rec.id, rec.taxID, rec.name));
			}
			records.clear();
		}

		outstream.println("\nTotal records: "+totalRecords+", condensed: "+condensed.size());
		outstream.println("Blacklisted bucket winners replaced: "+blacklisted);
		outstream.println("Both candidates blacklisted (kept higher): "+kept);

		validate(condensed);

		t.stop();
		outstream.println("\nTime: \t"+t);
	}

	private static LongHashSet getBlacklistSet(){
		try{
			java.lang.reflect.Field f=DynamicDemiLog.class.getDeclaredField("blacklist");
			f.setAccessible(true);
			return (LongHashSet)f.get(null);
		}catch(Exception e){
			throw new RuntimeException("Cannot access blacklist field", e);
		}
	}

	/*--------------------------------------------------------------*/
	/*-----------        Blacklist Generation        ---------------*/
	/*--------------------------------------------------------------*/

	void processBlacklist(Timer t){
		final TaxTree tree=TaxTree.sharedTree();
		final String[] inFiles=in.split(",");

		final LongObjectMap<IntHashSet> kmerGenera=new LongObjectMap<IntHashSet>();
		long totalRecords=0, totalWithKmers=0;
		int skippedTids=0;

		for(String inFile : inFiles){
			outstream.println("Loading "+inFile);
			ArrayList<DDLRecord> records=DDLLoader.loadFile(inFile.trim(), k);
			outstream.println("  "+records.size()+" records");

			int withKmers=0;
			for(DDLRecord rec : records){
				if(rec.ddl.hasKmers()){withKmers++;}
			}
			totalRecords+=records.size();
			totalWithKmers+=withKmers;

			for(DDLRecord rec : records){
				long[] kmers=rec.ddl.kmerArray();
				if(kmers==null){continue;}
				int tid=rec.taxID;
				if(tid<1 || (tree!=null && tid>=tree.nodes.length)){skippedTids++; continue;}
				int genusTid=tid;
				if(tree!=null){
					TaxNode gn=tree.getNodeAtLevel(tid, TaxTree.GENUS);
					if(gn!=null){genusTid=gn.id;}
				}
				for(long kmer : kmers){
					if(kmer==0){continue;}
					IntHashSet genera=kmerGenera.get(kmer);
					if(genera==null){
						genera=new IntHashSet(4);
						kmerGenera.put(kmer, genera);
					}
					genera.add(genusTid);
				}
			}
			records.clear();
			outstream.println("  Kmer map size: "+kmerGenera.size());
		}

		outstream.println("\nTotal records: "+totalRecords+" ("+totalWithKmers+" with kmers)");
		if(skippedTids>0){outstream.println("Skipped "+skippedTids+" records with out-of-range TaxIDs.");}
		outstream.println("Unique kmers across all files: "+kmerGenera.size());

		if(validate && inFiles.length==1){
			outstream.println("\nReloading for validation...");
			ArrayList<DDLRecord> records=DDLLoader.loadFile(inFiles[0].trim(), k);
			validate(records);
		}else if(validate){
			outstream.println("Validate mode requires a single input file; skipping.");
		}

		if(out!=null && totalWithKmers>0){
			buildBlacklist(kmerGenera, t);
		}else if(out!=null){
			outstream.println("ERROR: No records have kmer arrays. Rebuild with kmers=t.");
		}
	}

	/*--------------------------------------------------------------*/

	void buildBlacklist(LongObjectMap<IntHashSet> kmerGenera, Timer t){
		final TaxTree tree=TaxTree.sharedTree();
		final int[] LEVELS={TaxTree.GENUS, TaxTree.FAMILY, TaxTree.ORDER,
			TaxTree.CLASS, TaxTree.PHYLUM, TaxTree.KINGDOM, TaxTree.SUPERKINGDOM};
		final String[] LCODES={"g", "f", "o", "c", "p", "k", "sk"};

		final ArrayList<long[]> results=new ArrayList<long[]>();
		int byEntropy=0;
		final IntHashSet promoted=new IntHashSet(256);
		final long[] keys=kmerGenera.keys();
		final long invalid=kmerGenera.invalid();
		for(long kmer : keys){
			if(kmer==invalid){continue;}
			final IntHashSet genera=kmerGenera.get(kmer);
			final int genusCount=genera.size();
			final float ent=kmerEntropy(kmer, k, entropyK);
			final boolean passesTaxa=(genusCount>=minTaxCount);
			final boolean passesEntropy=(minEntropy>0 && ent<minEntropy);
			if(!passesTaxa && !passesEntropy){continue;}
			if(!passesTaxa){byEntropy++;}

			long[] entry=new long[2+LEVELS.length+1];
			entry[0]=kmer;
			entry[1]=genusCount;
			entry[2]=genusCount;
			entry[2+LEVELS.length]=(long)(ent*1000);

			if(tree!=null){
				int[] tids=genera.toArray();
				for(int lev=1; lev<LEVELS.length; lev++){
					promoted.clear();
					for(int tid : tids){
						TaxNode node=tree.getNodeAtLevel(tid, LEVELS[lev]);
						if(node!=null){promoted.add(node.id);}
					}
					entry[2+lev]=promoted.size();
				}
			}
			results.add(entry);
		}

		results.sort((a, b)->Long.compare(b[1], a[1]));
		outstream.println("Kmers output: "+results.size()
			+" ("+(results.size()-byEntropy)+" by genus count >= "+minTaxCount
			+(byEntropy>0 ? ", "+byEntropy+" by low entropy < "+minEntropy : "")
			+")");

		final FileFormat ff=FileFormat.testOutput(out, FileFormat.FASTA, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		for(long[] entry : results){
			final long kmer=entry[0];
			StringBuilder sb=new StringBuilder(80);
			sb.append(">kmer_").append(Long.toHexString(kmer));
			sb.append(" raw=").append(entry[1]);
			if(tree!=null){
				for(int lev=0; lev<LEVELS.length; lev++){
					sb.append(' ').append(LCODES[lev]).append('=').append(entry[2+lev]);
				}
			}
			sb.append(String.format(" e=%.3f", entry[2+LEVELS.length]/1000.0));
			bsw.print(sb.append('\n').toString());
			bsw.print(DynamicDemiLog.unpackKmer(kmer, k)+"\n");
		}
		bsw.poisonAndWait();

		t.stop();
		outstream.println("Wrote "+results.size()+" kmers to "+out);
		outstream.println("Time: \t"+t);
	}

	/*--------------------------------------------------------------*/

	void validate(ArrayList<DDLRecord> records){
		final TaxTree tree=TaxTree.sharedTree();
		if(tree==null){
			outstream.println("WARNING: TaxTree not loaded; cannot validate.");
			return;
		}

		final int[] levels={TaxTree.SPECIES, TaxTree.GENUS, TaxTree.FAMILY,
			TaxTree.ORDER, TaxTree.CLASS, TaxTree.PHYLUM, TaxTree.KINGDOM, TaxTree.SUPERKINGDOM};
		final String[] levelNames={"species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"};

		final shared.Random rand=Shared.threadLocalRandom(42);
		final int n=records.size();
		final long[] pairCounts=new long[levels.length];
		final double[] sharedSums=new double[levels.length];
		long skippedNoRank=0;

		outstream.println("\nValidating shared keys across taxonomic distances ("+validateSamples+" random pairs)...");
		for(int s=0; s<validateSamples; s++){
			int i=rand.nextInt(n), j=rand.nextInt(n);
			if(i==j){continue;}
			DDLRecord a=records.get(i), b=records.get(j);
			if(a.taxID<1 || b.taxID<1){continue;}
			if(a.taxID>=tree.nodes.length || b.taxID>=tree.nodes.length){continue;}

			TaxNode an=tree.getNode(a.taxID), bn=tree.getNode(b.taxID);
			if(an==null || bn==null){continue;}
			if(an.id==bn.id){continue;}
			TaxNode ca=tree.commonAncestor(an, bn);
			if(ca==null){continue;}
			final int caLevel=ca.level;

			int bucket=-1;
			for(int lev=0; lev<levels.length; lev++){
				if(caLevel==levels[lev]){bucket=lev; break;}
			}
			if(bucket<0){skippedNoRank++; continue;}
			if(bucket==0){continue;}

			int[] comp=a.ddl.compareToDetailed(b.ddl);
			int shared=comp[1];
			pairCounts[bucket]++;
			sharedSums[bucket]+=shared;
		}

		outstream.println("\nShared keys by common ancestor level:");
		outstream.println("Ancestor_level\tPairs\tAvg_shared_keys");
		for(int lev=0; lev<levels.length; lev++){
			if(pairCounts[lev]>0){
				outstream.println(levelNames[lev]+"\t"+pairCounts[lev]+"\t"+
					String.format("%.2f", sharedSums[lev]/pairCounts[lev]));
			}
		}
		if(skippedNoRank>0){
			outstream.println("(skipped "+skippedNoRank+" pairs with no-rank common ancestor)");
		}
		long unrelatedPairs=0;
		double unrelatedShared=0;
		for(int lev=levels.length/2; lev<levels.length; lev++){
			unrelatedPairs+=pairCounts[lev];
			unrelatedShared+=sharedSums[lev];
		}
		if(unrelatedPairs>0){
			outstream.println("(noise floor: "+String.format("%.2f", unrelatedShared/unrelatedPairs)+
				" avg shared keys at class+ distance)");
		}
	}

	/*--------------------------------------------------------------*/

	static float kmerEntropy(long kmer, int k, int ek){
		final int numSubkmers=k-ek+1;
		if(numSubkmers<1){return 1;}
		final int subkmerSpace=1<<(2*ek);
		final int subkmerMask=subkmerSpace-1;
		final int[] counts=new int[subkmerSpace];
		for(int i=0; i<numSubkmers; i++){
			int sub=(int)((kmer>>(2*i))&subkmerMask);
			counts[sub]++;
		}
		double entropy=0;
		final double logN=Math.log(numSubkmers);
		for(int c : counts){
			if(c>0){
				double p=(double)c/numSubkmers;
				entropy-=p*Math.log(p);
			}
		}
		return (float)(entropy/logN);
	}

	/*--------------------------------------------------------------*/

	private String in;
	private String out;
	private String blacklistFile;
	private int k=19;
	private int minTaxCount=20;
	private float minEntropy=0;
	private int entropyK=3;
	private boolean overwrite=false;
	private boolean verbose=false;
	private boolean validate=false;
	private boolean condense=false;
	private int validateSamples=100000;

	private static final PrintStream outstream=System.err;
}
