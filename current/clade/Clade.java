package clade;

import java.util.ArrayList;

import bin.AdjustEntropy;
import bin.SimilarityMeasures;
import prok.CallGenes;
import prok.GeneCaller;
import prok.Orf;
import shared.LineParser1;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import tax.PrintTaxonomy;
import tax.TaxNode;
import tax.TaxTree;
import tracker.EntropyTracker;
import tracker.KmerTracker;

/**
 * Represents a taxonomic clade with k-mer frequency signatures.
 * Contains 1-mer through 5-mer counts and derived statistics for genome comparison, using canonical k-mers to reduce dimensionality.
 * @author Brian Bushnell
 * @date April 12, 2025
 */
public class Clade extends CladeObject implements Comparable<Clade>{

	/**
	 * Constructs a Clade with the specified taxonomic information.
	 * Initializes k-mer count arrays from 1-mer through 5-mer using canonical indexing.
	 * @param taxID_ The taxonomic ID number
	 * @param level_ The taxonomic level (e.g., species, genus)
	 * @param name_ The taxonomic name
	 */
	public Clade(int taxID_, int level_, String name_) {
		taxID=taxID_;
		level=level_;
		name=name_;
		counts=new long[6][];
		counts[1]=new long[arrayLength[1]];
		counts[2]=new long[arrayLength[2]];
		counts[3]=new long[arrayLength[3]];
		counts[4]=new long[arrayLength[4]];
		counts[5]=new long[arrayLength[5]];//Could be optionally allocated later
	}

	public static Clade makeClade(int tid) {
		if(perSequence || !useTree) {return new Clade(tid, -1, null);}
		TaxNode tn=tree.getNode(tid);
		assert(tn!=null);
		if(tn==null) {
			return new Clade(tid, -1, null);//TODO: Changed to "Unknown_TaxID_" + tid in server version
		}
		return new Clade(tn.id, tn.level, tn.name);
	}
	
	public static Clade parseClade(ArrayList<byte[]> list, LineParser1 lp) {
		int lines=0;
		int coding=Clade.DECIMAL;
		int maxk=5;
		
		int pos=0;
		lp.set(list.get(pos));
		while(lp.startsWith('#')) {//Any number of header lines; but lines with tabs are parsed
			assert(lp.startsWith('#'));
			if(lp.terms()>1) {
				lines=lp.parseInt(1);
				for(int i=2; i<lp.terms(); i++) {
					if(lp.termEquals("DEC", i)) {coding=Clade.DECIMAL;}
					else if(lp.termEquals("A48", i)) {coding=Clade.A48;}
					else if(lp.termStartsWith("MAXK", i)) {maxk=lp.parseInt(i, 4);}
				}
			}
			pos++;
			lp.set(list.get(pos));
		}
		
		//First non-header line
		assert(lp.startsWith("tid"));
		int tid=lp.parseInt(1);
		
		pos++;
		lp.set(list.get(pos));
		assert(lp.startsWith("level"));
		int level=lp.parseInt(1);
		
		pos++;
		lp.set(list.get(pos));
		assert(lp.startsWith("name"));
		String name=lp.parseString(1);
		
		Clade c=new Clade(tid, level, name);
		
		pos++;
		synchronized(c) {
			if(Tools.startsWith(list.get(pos), "lineage")) {//Optional
				lp.set(list.get(pos));
				c.lineage=lp.parseString(1);
				pos++;
			}
			if(Tools.startsWith(list.get(pos), "gc")) {pos++;}//Optional
			if(Tools.startsWith(list.get(pos), "entropy")) {//Optional and slow, but should be in ref
				lp.set(list.get(pos));
				c.entropy=lp.parseFloat(1);
				pos++;
			}
			if(Tools.startsWith(list.get(pos), "strandedness")) {pos++;}//Calculated from dimers
			
			
			lp.set(list.get(pos));
			assert(lp.startsWith("bases")) : lp;//Could be calculated from monomers
			c.bases=lp.parseInt(1);

			pos++;
			lp.set(list.get(pos));
			assert(lp.startsWith("contigs")) : lp;
			c.contigs=lp.parseInt(1);

			for(int k=1; k<=maxk; k++) {
				//TODO: Change to concise packing, skipping k4/k5 depending on #bases, leaving them empty.
				pos++;
				lp.set(list.get(pos));
				assert(lp.startsWith((char)(k+'0')));
				final int terms=lp.terms(), expected=arrayLength[k]+1;
				assert(k<2 || terms==1 || terms==canonicalKmers[k]+1) : 
					k+", "+c.counts[k].length+", "+canonicalKmers[k]+", "+terms+", "+expected+", "+c.bases;
				assert(terms==expected || (k>3 && terms==1)) :
					k+", "+c.counts[k].length+", "+canonicalKmers[k]+", "+terms+", "+expected+", "+c.bases;
				
				if(terms>=expected) {
					if(coding==Clade.DECIMAL) {lp.parseLongArray(1, c.counts[k]);}
					else{lp.parseLongArrayA48(1, c.counts[k]);}
				}else{
					c.counts[k]=null;
				}
			}
			
			for(pos++; pos<list.size(); pos++) {
				lp.set(list.get(pos));
				if(lp.startsWith("16S")) {
					c.r16S=lp.parseByteArray(1);
				}else if(lp.startsWith("18S")) {
					c.r18S=lp.parseByteArray(1);
				}else if(lp.startsWith("k") && MAXK<5) {
					//do nothing
				}else {
					assert(false) : "Unknown line for TaxID "+c+"\n"+new String(list.get(pos));
				}
			}
			c.finish();
		}
		return c;
	}
	
	public synchronized void add(Read r, EntropyTracker et, GeneCaller caller) {
		add(r.bases, et);
		if(caller==null || hasSSU() || r.length()<900) {return;}
		assert(callSSU);
		ArrayList<Orf> genes=caller.callGenes(r);
		if(genes==null || genes.isEmpty()) {return;}
		for(Orf orf : genes) {
			if(orf.is16S()) {
				r16S=CallGenes.fetch(orf, r).bases;
				return;
			}else if(orf.is18S()) {
				r18S=CallGenes.fetch(orf, r).bases;
				return;
			}
		}
	}
	
	public synchronized void add(byte[] seq, EntropyTracker et) {
		finished=false;
		countKmersMulti(seq, counts, 5);
		float seqEntropy=(calcCladeEntropy ? et.averageEntropy(seq, false) : 0);
		entropy=(entropy*bases+seqEntropy*seq.length)/(float)(bases+seq.length);
		
		bases+=seq.length;
		contigs++;
	}
	
	/**
	 * Merges another Clade into this one, combining counts and updating statistics.
	 * If this Clade is empty, it adopts the taxonomic information of the other Clade.
	 * @param c The Clade to merge into this one
	 */
	public synchronized void add(Clade c) {
		finished=false;
		synchronized(c) {
			assert(c.taxID>0) : "\n"+this+"\n"+c+"\n";
			assert(c.bases>0) : "\n"+this+"\n"+c+"\n";//Not really necessary, but preventable
			assert(taxID==c.taxID || (taxID<1 && bases==0)) : "\n"+this+"\n"+c+"\n";
			if(taxID!=c.taxID) {
				assert(taxID<0 && bases==0);
				taxID=c.taxID;
				level=c.level;
				name=c.name;
				lineage=c.lineage;
			}
			
			Tools.add(counts, c.counts);
			entropy=(entropy*bases+c.entropy*c.bases)/(float)(bases+c.bases);

			bases+=c.bases;
			contigs+=c.contigs;
		}
	}
	
	/**
	 * Completes the Clade by calculating derived statistics.
	 * Computes GC content, strandedness, homopolymer and CA-GA measures, entropy compensation, and k-mer frequency arrays.
	 * Once completed, the Clade's state should not be modified.
	 */
	public synchronized void finish() {
		if(finished) {return;}
		gc=calcGC();
		strandedness=EntropyTracker.strandedness(counts[2], 2);
		hh=KmerTracker.HH(counts[2]);
		caga=KmerTracker.CAGA(counts[2]);
		gcCompEntropy=AdjustEntropy.compensate(gc, entropy);
		frequencies=new float[6][];
		frequencies[3]=toFrequencies(counts[3], 3);
		if(MAKE_FREQUENCIES && (method==ABSCOMP || method==ABS)) {
			frequencies[4]=(maxK<4 || bases<Comparison.minK4Bases ? null : toFrequencies(counts[4], 4));
			frequencies[5]=(maxK<5 || bases<Comparison.minK5Bases ? null : toFrequencies(counts[5], 5));
			if(DELETE_COUNTS) {counts[3]=counts[4]=counts[5]=null;}
		}
		finished=true;
	}
	
	private synchronized float calcGC() {
		long[] acgtn=counts[1];
		long a=acgtn[0], c=acgtn[1], g=acgtn[2], t=acgtn[3];
		return (float)((g+c)/Math.max(a+g+c+t, 1.0));
	}
	
	/**
	 * Makes a k-mer frequency array with raw or GC-compensated frequencies.
	 * For ABSCOMP, calls SimilarityMeasures.compensate to normalize within GC content groups; otherwise normalizes counts to sum to 1.
	 * Assumes counts come from a canonical k-mer array for compensated mode.
	 */
	public static synchronized float[] toFrequencies(long[] counts, final int k) {
		if(Comparison.method==Comparison.ABSCOMP) {return SimilarityMeasures.compensate(counts, k);}
		long sum=Tools.sum(counts);
		float inv=1f/sum;
		float[] freqs=new float[counts.length];
		for(int i=0; i<counts.length; i++) {freqs[i]=counts[i]*inv;}
		return freqs;
	}
	
	synchronized boolean hasSSU() {return r16S!=null | r18S!=null;}
	
	public synchronized void clear() {
		finished=false;
		
		taxID=level=-1;
		name=lineage=null;
		
		bases=contigs=0;
		gc=entropy=gcCompEntropy=strandedness=hh=caga=0;
		Tools.fill(counts, 0);
	}
	
	@Override
	public int compareTo(Clade b) {
		if(hh!=b.hh) {return hh>b.hh ? 1 : -1;}
		if(gc!=b.gc) {return gc>b.gc ? 1 : -1;}
		if(bases!=b.bases) {return bases>b.bases ? 1 : -1;}
		return taxID-b.taxID;
	}

	public CharSequence lineage() {
		if(lineage!=null) {return lineage;}
		return taxID<1 ? "NA" : (lineage=lineage(taxID).toString());
	}
	
	public static CharSequence lineage(int tid) {
		if(tree==null || tid<1) {return "NA";}
		TaxNode tn=tree.getNode(tid);
		if(tn==null) {return "NA";}
		return PrintTaxonomy.makeTaxLine(tree, tn, MIN_LINEAGE_LEVEL_E, TaxTree.SUPERKINGDOM_E, true, true);
	}
	
	public synchronized String toString() {
		ByteBuilder bb=new ByteBuilder();
		bb.append("tid=").append(taxID).append("\tgc=").append(gc, 4).append("\tname=").append(name);
		return bb.toString();
	}
	
	public synchronized ByteBuilder toBytes(ByteBuilder bb) {
		if(bb==null) {bb=new ByteBuilder();}
		assert(finished);
		if(!finished) {finish();}
		final boolean outDEC=(outputCoding==DECIMAL);
//		final byte[] temp=(outDEC ? null : KillSwitch.allocByte1D(12));
		{//header
			int lines=10+counts.length-1;
			if(r16S!=null) {lines++;}
			else if(r18S!=null) {lines++;}
			bb.append('#').tab().append(lines);
			bb.tab().append(outputCoding==DECIMAL ? "DEC" : "A48");
			if(MAXK!=5) {bb.tab().append("MAXK").append(MAXK);}
			bb.nl();
		}
		bb.append("tid\t").append(taxID).nl();
		bb.append("level\t").append(level).tab().append(TaxTree.levelToString(Math.max(0, level))).nl();
		bb.append("name\t").append(name).nl();
		if(writeLineage && taxID>1) {bb.append("lineage\t").append(lineage()).nl();}
		bb.append("gc\t").append(gc, 4).nl();
		bb.append("entropy\t").append(entropy, 8).nl();
		bb.append("strandedness\t").append(strandedness, 8).nl();
		bb.append("bases\t").append(bases).nl();
		bb.append("contigs\t").append(contigs).nl();
		int maxK=Math.min(MAXK, bases<Comparison.minK4Bases ? 3 : bases<Comparison.minK5Bases ? 4 : 5);
//		assert(false) : MAXK+", "+bases+", "+outDEC+", "+Arrays.toString(counts);
		for(int k=1; k<counts.length && k<=MAXK; k++) {
			bb.append(k).append("mers");
			if(counts[k]!=null && (k<=maxK || !CONCISE)) {
				bb.tab();
				if(outDEC) {bb.append(counts[k], '\t');}
//				else {bb.appendA48(counts[k], '\t', temp);}
				else {bb.appendA48(counts[k], '\t');}
			}
			bb.nl();
		}
		if(r16S!=null) {bb.append("16S\t").append(r16S).nl();}
		else if(r18S!=null) {bb.append("18S\t").append(r18S).nl();}
		return bb;
	}
	
	public synchronized boolean finished() {return finished;}
	
	public int taxID=-1;
	public int level=-1;
	public String name=null;
	public String lineage=null;
	public final long[][] counts;
	public float[][] frequencies;

	public byte[] r16S;
	public byte[] r18S;
	
	public long bases;
	public long contigs;
	public float gc;
	public float entropy;
	public float gcCompEntropy;
	public float strandedness;
	public float hh;
	public float caga;
	private boolean finished=false;
	
	public static final int DECIMAL=0, A48=1;
	public static int outputCoding=A48; //A48 breaks Cloudflare...  now bypassed
	public static int MAXK=5;
	public static boolean callSSU=false;
	public static boolean writeLineage=true;
	public static boolean MAKE_FREQUENCIES=true;
	public static boolean DELETE_COUNTS=false;//Only OK when searching local index. Which includes on server.
	public static boolean CONCISE=true;//TODO: Set to true once tested and running on server
	
}