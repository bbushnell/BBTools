package dna;

import java.util.ArrayList;

import shared.Shared;
import shared.Tools;

/**
 * A group of transcript {@link Gene}s intended to share one gene id and chromosome.
 * The caller supplies the group name; it is not derived from each transcript.
 * Aggregates their span (minStart..maxEnd) and folds per-transcript flags: pseudo and
 * untranslated are the logical AND across all transcripts, and strand collapses to '?'
 * if the transcripts disagree. Ordered/compared by chromosome, then minStart, then name.
 * The transcript list is borrowed and public; summaries are computed only at construction
 * and are not refreshed after list edits. Mutable span fields must remain unchanged
 * while an instance participates in an ordered collection. No synchronization is provided.
 * @author Brian Bushnell
 */
public class GeneSet implements Comparable<GeneSet>{

	/** Loads the configured gene-ID table as a legacy diagnostic; ignores arguments. */
	public static void main(String[] args){
		Data.getGeneIDTable();
	}

	/** Retains a nonempty list of nonnull transcripts and computes its aggregate fields.
	 * All transcript IDs must agree; chromosome is taken from the first transcript.
	 * Requires a nonnull name for comparison/hashing and valid transcript coordinates. */
	public GeneSet(String n, ArrayList<Gene> g){
		name=n;
		id=g.get(0).id;
		genes=g;
		chrom=g.get(0).chromosome;
		transcripts=genes.size();
		assert(transcripts>0);

		byte st=-1;

		boolean pse=true, unt=true;

		for(int i=0; i<transcripts; i++){
			Gene gene=g.get(i);

			assert(gene.id==id) : g;

			pse=(pse && gene.pseudo);
			unt=(unt && gene.untranslated);
			minStart=min((int)gene.txStart, minStart);
			maxEnd=max((int)gene.txStop, maxEnd);
			//				assert(st==-1 || st==gene.strand) : g;
			if(st==-1){st=gene.strand;}
			else if(st!=gene.strand){st=(byte) Tools.find("?", Shared.strandCodes);}
		}

		pseudo=pse;
		untranslated=unt;

		for(Gene gene : g){
			assert(pseudo==gene.pseudo || (!pseudo && gene.pseudo)) : g;
			assert(untranslated==gene.untranslated || (!untranslated && gene.untranslated)) : g;
//			assert(untranslated==gene.untranslated) : g;
		}

		strand=st;
	}

	/** Formats the name, transcript count, inclusive aggregate span and strand symbol. */
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append(name);
		while(sb.length()<10){sb.append(' ');}
		sb.append('\t');
		sb.append(padFront(transcripts+"", 2)+" transcript"+(transcripts==1 ? " " : "s"));

		sb.append("\tchr"+chrom+" ("+minStart+" - "+maxEnd+"), '"+Shared.strandCodes[strand]+"'");

		return sb.toString();
	}

	/** Left-pads with spaces to at least width characters; never truncates. */
	private static final String padFront(String num, int width){
		String r=num;
		while(r.length()<width){r=" "+r;}
		return r;
	}

	/** Right-pads with spaces to at least width characters; never truncates. */
	private static final String padBack(String num, int width){
		String r=num;
		while(r.length()<width){r=r+" ";}
		return r;
	}

	/** Caller-supplied grouping name used in comparison and hashing. */
	public final String name;
	/** ID of the first transcript; constructor assertions require all IDs to match. */
	public final int id;
	/** Chromosome of the first transcript. */
	public final byte chrom;
	/** Common strand, or the '?' code when valid transcript strands disagree. */
	public final byte strand;
	/** Borrowed mutable transcript list; changes do not refresh aggregate fields. */
	public final ArrayList<Gene> genes;
	/** List size captured at construction. */
	public final int transcripts;

	/** True only if every transcript is untranslated at construction. */
	public final boolean untranslated;
	/** True only if every transcript is a pseudogene at construction. */
	public final boolean pseudo;

	/** Minimum inclusive transcription start, initially reduced from Integer.MAX_VALUE. */
	public int minStart=Integer.MAX_VALUE;
	/** Maximum inclusive transcription stop, initially reduced from zero. */
	public int maxEnd=0;

	/** Whether a coordinate lies inside the aggregate span, including intronic gaps. */
	public boolean intersects(int point){
		return point>=minStart && point<=maxEnd;
	}
	/** Whether an ordered inclusive interval overlaps the aggregate span. */
	public boolean intersects(int point1, int point2){
		return point2>=minStart && point1<=maxEnd;
	}

	/** Orders by chromosome, minimum start and name; ignores transcript contents and end. */
	@Override
	public int compareTo(GeneSet other){
		if(chrom!=other.chrom){
			return chrom>other.chrom ? 1 : -1;
		}
		int x=minStart<other.minStart ? -1 : minStart>other.minStart ? 1 : 0;
		if(x!=0){return x;}
		return x=name.compareTo(other.name);
	}

	/**
	 * Tests equality by casting to GeneSet and delegating to {@code equals(GeneSet)}.
	 * Expects {@code other} to be a GeneSet and will throw {@link ClassCastException} otherwise.
	 * The current typed implementation also throws for null.
	 * @param other Object to compare (must be a GeneSet)
	 * @return true if objects are equal GeneSets
	 */
	@Override
	public boolean equals(Object other){
		return equals((GeneSet)other);
	}

	/** Equality of chromosome, minimum start and name; currently requires a nonnull argument. */
	public boolean equals(GeneSet other){
		//Note: equality is by compareTo, i.e. chromosome+minStart+name only — NOT gene content.
		//Two GeneSets with the same name/chrom/minStart but different transcripts compare equal.
		//Consistent with hashCode (name+chrom).
		return compareTo(other)==0;
	}

	/** Hashes name and chromosome; unequal starts may legitimately share the hash. */
	@Override
	public int hashCode(){
		return Integer.rotateLeft(name.hashCode(), 5)^chrom;
	}

	/** Smaller integer. */
	private static final int min(int x, int y){return x<y ? x : y;}
	/** Larger integer. */
	private static final int max(int x, int y){return x>y ? x : y;}
	/** Smaller byte. */
	private static final byte min(byte x, byte y){return x<y ? x : y;}
	/** Larger byte. */
	private static final byte max(byte x, byte y){return x>y ? x : y;}
	/** Smaller long. */
	private static final long min(long x, long y){return x<y ? x : y;}
	/** Larger long. */
	private static final long max(long x, long y){return x>y ? x : y;}
	/** Returns x when x is less than y, otherwise y; does not special-case NaN. */
	private static final float min(float x, float y){return x<y ? x : y;}
	/** Returns x when x is greater than y, otherwise y; does not special-case NaN. */
	private static final float max(float x, float y){return x>y ? x : y;}

}
