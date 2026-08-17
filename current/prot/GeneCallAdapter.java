package prot;

import java.util.ArrayList;

import dna.Data;
import fileIO.FileFormat;
import prok.CallGenes;
import prok.GeneCaller;
import prok.GeneModel;
import prok.Orf;
import prok.PGMTools;
import prok.ProkObject;
import stream.Read;
import stream.ReadInputStream;

/**
 * In-memory bridge from raw contigs to predicted proteins for the MAG-QC pipeline.
 *
 * <p>Drives the existing prokaryotic gene caller ({@link prok.GeneCaller}) on contigs
 * already resident in memory and returns the translated CDS as a list of
 * {@link ProteinSequence}, with no intermediate disk file. This is the same
 * gene-calling and translation path used by {@code callgenes.sh outa=}, so results
 * match that tool's amino-acid FASTA; the only new capability is that it stays
 * in memory (contigs in, {@code List<ProteinSequence>} out).</p>
 *
 * <p>Gene calling is scoped to CDS only. The rRNA/tRNA callers require consensus
 * sequences and long-kmer sets that the CallGenes driver loads separately; the
 * MAG-QC protein path does not need them, and disabling them also avoids running
 * the RNA callers with unloaded reference data. The relevant ProkObject call-flags
 * are saved and restored around each call so global state is left unchanged.</p>
 *
 * <p><b>Thread safety:</b> {@link ProkObject}'s call-flags are shared mutable
 * statics, read repeatedly throughout {@link GeneCaller#callGenes}'s per-contig,
 * per-frame execution (not just once at entry) -- so two threads calling
 * {@link #callProteins(ArrayList,GeneModel)} concurrently could interleave their
 * save/set/restore sequences, with one call's restore firing while another is
 * still mid-execution and expecting the CDS-only scoping to hold. Demonstrated
 * 2026-08-16: 30 rounds of one long call racing many short calls left the global
 * RNA-call flags permanently stuck (not restored to baseline) in 7/30 rounds.
 * {@link #callProteins(ArrayList,GeneModel)} is therefore {@code synchronized} --
 * the save/set/scoped-work/restore sequence is atomic with respect to any other
 * concurrent caller of this class. This serializes concurrent callers of this
 * adapter (a real cost if it's ever driven from a per-bin thread pool), but the
 * underlying flags are {@link ProkObject} statics shared far beyond this class,
 * so making them thread-local instead would be a larger, separately-scoped
 * change to core prok/ infrastructure -- flagged, not done here.</p>
 *
 * @author Eru
 */
public final class GeneCallAdapter {

	/** Utility class; use the static methods. */
	private GeneCallAdapter(){}

	/** Cached default gene model (resources/model.pgm), loaded on first use. */
	private static GeneModel defaultModel;

	/**
	 * Loads and caches the shipped default gene model (resources/model.pgm).
	 * @return The default {@link GeneModel}, loaded via {@link Data#findPath} the first time.
	 */
	public static synchronized GeneModel defaultModel(){
		if(defaultModel==null){
			final String path=Data.findPath("?model.pgm");
			assert(path!=null) : "Could not locate the default model.pgm via Data.findPath(\"?model.pgm\").";
			final ArrayList<String> list=new ArrayList<String>(1);
			list.add(path);
			defaultModel=PGMTools.loadAndMerge(list);
			assert(defaultModel!=null) : "PGMTools.loadAndMerge returned null for "+path;
		}
		return defaultModel;
	}

	/**
	 * Calls CDS and returns their translated proteins, using the shipped default model.
	 * @param contigs Genome/bin contigs already in memory (nucleotide {@link Read}s).
	 * @return Predicted proteins as {@link ProteinSequence}; empty if no genes are called.
	 */
	public static ArrayList<ProteinSequence> callProteins(final ArrayList<Read> contigs){
		return callProteins(contigs, defaultModel());
	}

	/**
	 * Calls CDS and returns their translated proteins, using the supplied model.
	 *
	 * <p>Each contig is passed through {@link GeneCaller#callGenes(Read, GeneModel, boolean)}
	 * and the resulting CDS ORFs are translated with {@link CallGenes#translate(Read, ArrayList)}
	 * (the terminal stop codon is dropped by the translator). Each protein's id is the
	 * source contig id plus strand and coordinates, tab-normalized to be a valid
	 * {@link ProteinSequence} identifier.</p>
	 *
	 * @param contigs Genome/bin contigs already in memory (nucleotide {@link Read}s).
	 * @param pgm Gene model to score ORFs against.
	 * @return Predicted proteins as {@link ProteinSequence}; empty if no genes are called.
	 */
	public static synchronized ArrayList<ProteinSequence> callProteins(final ArrayList<Read> contigs, final GeneModel pgm){
		assert(contigs!=null);
		assert(pgm!=null);

		//Save global RNA call-flags, then scope this call to CDS only.
		final boolean oldTRNA=ProkObject.calltRNA, old16S=ProkObject.call16S, old23S=ProkObject.call23S;
		final boolean old5S=ProkObject.call5S, old18S=ProkObject.call18S, oldCDS=ProkObject.callCDS;
		ProkObject.callCDS=true;
		ProkObject.calltRNA=ProkObject.call16S=ProkObject.call23S=ProkObject.call5S=ProkObject.call18S=false;

		final ArrayList<ProteinSequence> out=new ArrayList<ProteinSequence>();
		try{
			final GeneCaller caller=CallGenes.makeGeneCaller(pgm);
			for(final Read r : contigs){
				if(r==null || r.bases==null){continue;}
				final ArrayList<Orf> orfs=caller.callGenes(r, pgm, true);
				final ArrayList<Read> prots=CallGenes.translate(r, orfs);
				if(prots==null){continue;}
				for(final Read p : prots){
					if(p==null || p.bases==null || p.bases.length==0){continue;}
					final String id=sanitizeId(p.id);
					out.add(new ProteinSequence(id, p.bases));
				}
			}
		}finally{
			//Restore global call-flags no matter what.
			ProkObject.callCDS=oldCDS;
			ProkObject.calltRNA=oldTRNA;
			ProkObject.call16S=old16S;
			ProkObject.call23S=old23S;
			ProkObject.call5S=old5S;
			ProkObject.call18S=old18S;
		}
		return out;
	}

	/**
	 * Convenience: loads a genome/bin FASTA from disk into memory and calls proteins.
	 * Provided for the CLI and testing; the in-memory {@link #callProteins(ArrayList)}
	 * overloads are the real new capability.
	 * @param fastaPath Path to a nucleotide FASTA file (may be gzipped).
	 * @return Predicted proteins as {@link ProteinSequence}.
	 */
	public static ArrayList<ProteinSequence> callProteins(final String fastaPath){
		assert(fastaPath!=null);
		final ArrayList<Read> contigs=ReadInputStream.toReads(fastaPath, FileFormat.FASTA, -1);
		assert(contigs!=null && !contigs.isEmpty()) : "No contigs read from "+fastaPath;
		return callProteins(contigs);
	}

	/**
	 * Normalizes a translated-protein Read id into a valid ProteinSequence identifier.
	 * The translator encodes coordinates as tab-delimited fields ("id\tstrand\tstart-stop");
	 * tabs are illegal in a {@link ProteinSequence} id, so they are replaced with '_'.
	 * @param rawId The Read id produced by the translator.
	 * @return A tab-free identifier.
	 */
	private static String sanitizeId(final String rawId){
		assert(rawId!=null && rawId.length()>0);
		return rawId.replace('\t', '_');
	}

	/**
	 * CLI: genome FASTA in, protein summary (and optional protein FASTA) out.
	 * Usage: {@code java prot.GeneCallAdapter <genome.fasta> [out.faa]}
	 * Note: {@code callgenes.sh in=genome outa=out.faa} already covers file-to-file
	 * protein output; this CLI exists mainly to exercise the in-memory adapter.
	 * @param args Command-line arguments: input FASTA, optional output FASTA path.
	 */
	public static void main(final String[] args){
		if(args.length<1){
			System.err.println("Usage: java prot.GeneCallAdapter <genome.fasta> [out.faa]");
			System.exit(1);
		}
		final String in=args[0];
		final ArrayList<ProteinSequence> proteins=callProteins(in);

		System.err.println("Input:    "+in);
		System.err.println("Proteins: "+proteins.size());
		if(!proteins.isEmpty()){
			final ProteinSequence first=proteins.get(0);
			System.err.println("First id: "+first.id);
			System.err.println("First len:"+first.length());
			final int show=Math.min(30, first.enc.length);
			final StringBuilder sb=new StringBuilder(show);
			for(int i=0; i<show; i++){
				sb.append(decode(first.enc[i]));
			}
			System.err.println("First 30: "+sb);
		}

		if(args.length>1){
			final String outPath=args[1];
			final fileIO.ByteStreamWriter bsw=new fileIO.ByteStreamWriter(outPath, true, false, true, FileFormat.FA);
			bsw.start();
			final structures.ByteBuilder bb=new structures.ByteBuilder();
			for(final ProteinSequence ps : proteins){
				bb.clear();
				bb.append('>').append(ps.id).nl();
				for(int i=0; i<ps.enc.length; i++){bb.append(decode(ps.enc[i]));}
				bb.nl();
				bsw.print(bb);
			}
			bsw.poisonAndWait();
			System.err.println("Wrote:    "+outPath);
		}
	}

	/**
	 * Decodes one Blosum62-encoded residue back to its one-letter ASCII code.
	 * Standard residues (0-19) invert through {@link dna.AminoAcid#numberToAcid};
	 * {@link Blosum62#X_CODE} maps to 'X'.
	 * @param enc Encoded residue value (0-19, or {@link Blosum62#X_CODE}).
	 * @return The one-letter amino-acid character.
	 */
	private static char decode(final byte enc){
		if(enc>=0 && enc<20){return (char)dna.AminoAcid.numberToAcid[enc];}
		return 'X';
	}
}
