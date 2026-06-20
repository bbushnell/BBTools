package bin;

import dna.Data;
import prok.CallGenes;
import prok.GeneCaller;
import prok.GeneModel;
import prok.GeneModelParser;

/**
 * Static utility class for configuring and managing gene detection models in prokaryotic and eukaryotic genomes.
 * Provides synchronized methods to initialize, configure, and retrieve gene callers using pre-defined Markov models.
 * Supports lazy initialization of gene detection models with automatic file path resolution.
 *
 * @author Brian Bushnell
 * @date December 5, 2024
 */
public class GeneTools {
	//=== Eru review 2026-06-20 (V3) === static gene-caller config/loader (lazy PGM load); feeds the QuickBin gene-calling path.
	//#001 LOW (latent): loadPGM:44 ternary `gCaller=(pgm==null && gCaller!=null ? null : CallGenes.makeGeneCaller(pgm))` is
	//  muddled — when pgm==null && gCaller==null (first call + model-load FAILURE) it calls makeGeneCaller(NULL) instead of
	//  leaving gCaller null; the `&& gCaller!=null` clause is meaningless. Intended: gCaller=(pgm==null?null:makeGeneCaller(pgm)).
	//  Latent (model.pgm ships; only bites if missing/corrupt). Not patched — error-path behavior change; verify makeGeneCaller(null)+intent w/ Brian.
	//NOTE: loadPGM sets CallGenes.call16S/call18S (40) but setMode sets GeneCaller.call16S/etc — same statics or different flag
	//  targets? possible mismatch (comprehension QUESTION, cross-package into prok/). pgmFile is set at static-init (70) so the
	//  `if(pgmFile==null)` re-resolve at loadPGM:37 is dead. All methods synchronized (thread-safe lazy init).

	/**
	 * Creates a new GeneCaller instance using the loaded gene model.
	 * Automatically loads the prokaryotic gene model if not already loaded.
	 * Thread-safe through synchronization.
	 * @return A new GeneCaller configured with the current gene model
	 */
	public static synchronized GeneCaller makeGeneCaller() {
		if(pgm==null) {loadPGM();}
		return CallGenes.makeGeneCaller(pgm);
	}
	
	/**
	 * Loads the prokaryotic gene model (PGM) from disk if not already loaded.
	 * Initializes long k-mers, consensus sequences, and enables 16S/18S detection.
	 * Thread-safe through synchronization with early return if model already loaded.
	 */
	public static synchronized void loadPGM() {
		if(pgm!=null) {return;}
		if(pgmFile==null){pgmFile=Data.findPath("?model.pgm");}
		
		if(!quiet) {System.err.println("Loading "+pgmFile);}
		CallGenes.call16S=CallGenes.call18S=true;
		CallGenes.loadLongKmers();
		CallGenes.loadConsensusSequenceFromFile(true, true);
		pgm=GeneModelParser.loadModel(pgmFile);
		gCaller=(pgm==null && gCaller!=null ? null : CallGenes.makeGeneCaller(pgm));
	}
	
	/**
	 * Configures gene detection modes for different RNA types and coding sequences.
	 * Sets which gene types should be detected during analysis.
	 *
	 * @param r16 Enable 16S ribosomal RNA detection
	 * @param r18 Enable 18S ribosomal RNA detection
	 * @param r5 Enable 5S ribosomal RNA detection
	 * @param r23 Enable 23S ribosomal RNA detection
	 * @param trna Enable transfer RNA detection
	 * @param cds Enable coding sequence detection
	 */
	static synchronized void setMode(boolean r16, boolean r18, boolean r5, boolean r23, boolean trna, boolean cds) {
		GeneCaller.call16S=r16;
		GeneCaller.call18S=r18;
		GeneCaller.call5S=r5;
		GeneCaller.call23S=r23;
		GeneCaller.calltRNA=trna;
		GeneCaller.callCDS=cds;
	}
	
	/**
	 * File path to the prokaryotic gene model, resolved automatically using Data.findPath
	 */
	static String pgmFile=Data.findPath("?model.pgm");
	/** Loaded prokaryotic gene model used for gene detection and classification */
	static GeneModel pgm;
	/** Cached GeneCaller instance created from the loaded gene model */
	static GeneCaller gCaller;
	/** Controls whether model loading progress messages are printed to stderr */
	static boolean quiet=false;
	
}
