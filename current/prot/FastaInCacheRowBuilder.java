package prot;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;

import dna.Data;
import prok.CallGenes;
import prok.GeneCaller;
import prok.GeneModel;
import prok.Orf;
import prok.ProkObject;
import stream.Read;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntList;
import tax.TaxTree;
import tracker.KmerTracker;

/**
 * Builds CacheBuilder's 19-field rows directly from in-memory nucleotide contigs.
 *
 * <p>This is the F2 production seam: the same CallGenes path is used with all
 * RNA callers enabled, feature statistics are accumulated through
 * {@link CacheBuilder.Acc#addOrf(Orf)}, and the final row is formatted by
 * {@link CacheBuilder#appendRow}.  The optional family-assignment overload
 * translates those same CDS ORFs and calls {@link ProteinSearcher} directly;
 * no GFF, protein FASTA, or hits file is created.</p>
 *
 * @author Sayu
 */
public final class FastaInCacheRowBuilder {

	/** Utility class; use {@link #buildRows}. */
	private FastaInCacheRowBuilder(){}

	/**
	 * Builds one row per contig using the shipped default gene model.
	 *
	 * @param contigs Nucleotide contigs belonging to one domain.
	 * @param domain Either {@code Bacteria} or {@code Archaea}.
	 * @return CacheBuilder-formatted data rows, without the header or newline.
	 */
	public static synchronized ArrayList<String> buildRows(final ArrayList<Read> contigs,
			final String domain){
		requireLegacyCaller();
		return buildRows(contigs, domain, GeneCallAdapter.defaultModel(), 1);
	}

	/**
	 * Builds rows with the caller-supplied per-phylum gene model.  The model is
	 * required by the real training path; the no-model overload is only a
	 * convenience for callers using the shipped default model.
	 *
	 * @param contigs Nucleotide contigs belonging to one domain.
	 * @param domain Either {@code Bacteria} or {@code Archaea}.
	 * @param pgm Per-phylum (or explicitly chosen) gene model.
	 * @return CacheBuilder-formatted data rows, without the header or newline.
	 */
	public static synchronized ArrayList<String> buildRows(final ArrayList<Read> contigs,
			final String domain, final GeneModel pgm){
		return buildRows(contigs, domain, pgm, 1);
	}

	/** Builds rows with an explicit CallGenes self-refinement pass count. */
	public static synchronized ArrayList<String> buildRows(final ArrayList<Read> contigs,
			final String domain, final GeneModel pgm, final int passes){
		return buildRows(contigs, domain, pgm, passes, true);
	}

	/**
	 * Builds rows and assigns each translated CDS to ZERO-OR-ONE family via the canonical
	 * {@link ProteinSearcher#assignFamily} (Patch 2). No search()/considerHit reduction and NO
	 * fallback: a rejected CDS records no family copy. One family copy per ASSIGNED CDS is recorded
	 * at its roster rank ({@link ProteinSearcher.FamilyAssignment#familyIdx}) -- the same rank space
	 * the removed legacy path recorded via {@code state.rank}, so cache roster indexing is preserved.
	 * The binding, searcher (its {@code aligner} must be {@code "d55"}, assignFamily's own guard), and
	 * policy are all required. THE ONLY family-assignment entrypoint.
	 */
	public static synchronized ArrayList<String> buildRowsWithAssignments(final ArrayList<Read> contigs,
			final String domain, final GeneModel pgm, final int passes,
			final ProteinSearcher.AssignmentBinding binding, final ProteinSearcher searcher,
			final ProteinSearcher.AssignPolicy policy){
		return buildRowsWithAssignments(contigs,domain,pgm,passes,binding,searcher,policy,0);
	}

	/**
	 * Builds rows using an explicit bounded-lookahead distance. Other assignment
	 * policies require {@code lookahead==0}; invalid combinations fail before
	 * gene calling begins.
	 */
	public static synchronized ArrayList<String> buildRowsWithAssignments(final ArrayList<Read> contigs,
			final String domain, final GeneModel pgm, final int passes,
			final ProteinSearcher.AssignmentBinding binding, final ProteinSearcher searcher,
			final ProteinSearcher.AssignPolicy policy, final int lookahead){
		if(binding==null || searcher==null || policy==null){
			throw new RuntimeException("AssignmentBinding, searcher, and policy are all required.");
		}
		ProteinSearcher.validateLookahead(policy,lookahead);
		return buildRows(contigs,domain,pgm,passes,true,searcher,binding,policy,lookahead);
	}

	/**
	 * Builds one typed whole-bin observation without crossing the historical cache-row boundary.
	 * The caller supplies the already-resolved QuickClade fields and the exact gene model; this
	 * method never parses a TID, serializes a cache row, or derives taxonomy.  Family width is
	 * taken from the same immutable binding used by the canonical assignment API.
	 */
	static synchronized MagQCPreparedBin buildPreparedBin(final ArrayList<Read> contigs,
			final String binId, final String qcStatus, final String qcDomain, final String qcPhylum,
			final GeneModel pgm, final int passes, final ProteinSearcher.AssignmentBinding binding,
			final ProteinSearcher searcher, final ProteinSearcher.AssignPolicy policy, final int lookahead){
		return buildPreparedBin(contigs, binId, qcStatus, qcDomain, qcPhylum, pgm, passes,
			binding, searcher, policy, lookahead, false);
	}

	/** Shared arithmetic; only the caller lifecycle differs between legacy and batch operation. */
	private static MagQCPreparedBin buildPreparedBin(final ArrayList<Read> contigs,
			final String binId, final String qcStatus, final String qcDomain, final String qcPhylum,
			final GeneModel pgm, final int passes, final ProteinSearcher.AssignmentBinding binding,
			final ProteinSearcher searcher, final ProteinSearcher.AssignPolicy policy, final int lookahead,
			final boolean configured){
		if(binding==null || searcher==null || policy==null){
			throw new RuntimeException("Prepared assignment requires binding, searcher, and policy.");
		}
		ProteinSearcher.validateLookahead(policy, lookahead);
		if(contigs==null || contigs.isEmpty()){throw new RuntimeException("No contigs supplied.");}
		if(pgm==null){throw new RuntimeException("No gene model supplied.");}
		if(passes<1){throw new RuntimeException("CallGenes passes must be >=1: "+passes);}
		final CallerData data=configured ? callConfiguredContigs(contigs, pgm, passes, true) :
			callContigs(contigs, pgm, passes, true, true);
		assignFamiliesViaBinding(data.accs, data.proteins, data.proteinToContig,
			binding, searcher, policy, lookahead);
		final MagQCPreparedBin bin=new MagQCPreparedBin(binding.familyCount());
		bin.id=binId; bin.status=qcStatus; bin.domain=qcDomain; bin.phylum=qcPhylum;
		for(final Read r : contigs){
			if(r==null || r.bases==null){continue;}
			for(final byte ch : r.bases){
				bin.length=addPositive(bin.length, 1L, "length_bp", binId);
				switch(ch){
					case 'G': case 'g': case 'C': case 'c':
						bin.agg.gc=addPositive(bin.agg.gc, 1L, "gc_bases", binId);
						bin.agg.acgt=addPositive(bin.agg.acgt, 1L, "acgt_bases", binId); break;
					case 'A': case 'a': case 'T': case 't':
						bin.agg.acgt=addPositive(bin.agg.acgt, 1L, "acgt_bases", binId); break;
					default: break;
				}
			}
		}
		final KmerTracker binDimers=new KmerTracker(2);
		for(final Read r : contigs){
			if(r==null || r.bases==null){continue;}
			binDimers.add(r.bases);
		}
		for(int i=0; i<bin.agg.dimer.length; i++){
			final long total=binDimers.counts[i];
			if(total>Integer.MAX_VALUE){throw overflow("dimer["+i+"]", binId);}
			bin.agg.dimer[i]=(int)total;
		}
		for(final CacheBuilder.Acc acc : data.accs.values()){
			bin.agg.cds=addInt(bin.agg.cds, acc.cds, "cds", binId);
			bin.agg.mapped=addInt(bin.agg.mapped, acc.mapped, "mapped_cds", binId);
			bin.agg.coding=addPositive(bin.agg.coding, acc.coding, "coding_bp", binId);
			bin.agg.glenSum=addPositive(bin.agg.glenSum, acc.glenSum, "gene_length_sum", binId);
			bin.agg.glenSq=addPositive(bin.agg.glenSq, acc.glenSq, "gene_length_sq_sum", binId);
			bin.agg.r16=addInt(bin.agg.r16, acc.r16, "r16", binId);
			bin.agg.r23=addInt(bin.agg.r23, acc.r23, "r23", binId);
			bin.agg.r5=addInt(bin.agg.r5, acc.r5, "r5", binId);
			bin.agg.rother=addInt(bin.agg.rother, acc.rother, "rother", binId);
			bin.agg.trna=addInt(bin.agg.trna, acc.trna, "trna_total", binId);
			for(final int rank : acc.fam.keys()){
				if(rank==acc.fam.invalid()){continue;}
				final int count=acc.fam.get(rank);
				if(rank<0 || rank>=bin.families.length){throw new RuntimeException("Family rank out of prepared range: "+rank);}
				bin.families[rank]=addInt(bin.families[rank], count, "family_counts", binId);
			}
			final int invalid=acc.anticodon.invalid();
			for(final int code : acc.anticodon.keys()){
				if(code==invalid){continue;}
				final int count=acc.anticodon.get(code), old=bin.agg.anti[code];
				bin.agg.anti[code]=addInt(old, count, "anticodon_counts", binId);
			}
		}
		bin.validate();
		return bin;
	}

	/**
	 * Builds rows with an explicit boundary-refinement setting.  The
	 * package-visible switch exists for the focused contig-71 regression; normal
	 * production callers use the default-on four-argument overload.
	 */
	static synchronized ArrayList<String> buildRows(final ArrayList<Read> contigs,
			final String domain, final GeneModel pgm, final int passes, final boolean boundaryNet){
		return buildRows(contigs,domain,pgm,passes,boundaryNet,null,null,null,0);
	}

	private static ArrayList<String> buildRows(final ArrayList<Read> contigs,
			final String domain, final GeneModel pgm, final int passes, final boolean boundaryNet,
			final ProteinSearcher searcher, final ProteinSearcher.AssignmentBinding binding,
			final ProteinSearcher.AssignPolicy policy, final int lookahead){
		if(contigs==null || contigs.isEmpty()){throw new RuntimeException("No contigs supplied.");}
		if(!"Bacteria".equals(domain) && !"Archaea".equals(domain)){
			throw new RuntimeException("Unsupported cache domain: "+domain);
		}
		if(pgm==null){throw new RuntimeException("No gene model supplied.");}
		if(passes<1){throw new RuntimeException("CallGenes passes must be >=1: "+passes);}
		final CallerData data=callContigs(contigs, pgm, passes, boundaryNet, binding!=null);
		final HashMap<String, CacheBuilder.Acc> accs=data.accs;
		if(binding!=null){assignFamiliesViaBinding(accs,data.proteins,data.proteinToContig,
				binding,searcher,policy,lookahead);}

		final IntHashMap archaeaSet=new IntHashMap(Math.max(4, contigs.size()));
		final IntHashMap bacteriaSet=new IntHashMap(Math.max(4, contigs.size()));
		for(final String name : accs.keySet()){
			final int tid=TaxTree.parseTaxID(name);
			if(tid<=0){throw new RuntimeException("Could not parse positive tid from contig: "+name);}
			if("Archaea".equals(domain)){archaeaSet.put(tid, 1);}
			else{bacteriaSet.put(tid, 1);}
		}

		final ArrayList<String> rows=new ArrayList<String>(accs.size());
		final IntList rankBuf=new IntList(64), countBuf=new IntList(64);
		final KmerTracker dimerTracker=new KmerTracker(2);
		for(final Read r : contigs){
			if(r==null || r.bases==null){continue;}
			final String name=canonicalName(r.id);
			int length=0, gc=0, acgt=0;
			dimerTracker.clearAll();
			for(final byte ch : r.bases){
				length++;
				switch(ch){
					case 'G': case 'g': case 'C': case 'c': gc++; acgt++; break;
					case 'A': case 'a': case 'T': case 't': acgt++; break;
					default: break;
				}
				dimerTracker.add(ch);
			}
			final ByteBuilder bb=new ByteBuilder(256);
			CacheBuilder.appendRow(bb, accs, archaeaSet, bacteriaSet, name, length, gc, acgt,
				dimerTracker.counts, rankBuf, countBuf);
			final String row=bb.toString();
			rows.add(row.substring(0, row.length()-1));
		}
		return rows;
	}

	/** Shared CallGenes preparation and ORF extraction for cache rows and prepared bins. */
	private static CallerData callContigs(final ArrayList<Read> contigs, final GeneModel pgm,
			final int passes, final boolean boundaryNet, final boolean collectProteins){
		requireLegacyCaller();
		try(CallerSettings previous=new CallerSettings()){
			configureCaller(boundaryNet);
			return callConfiguredContigs(contigs, pgm, passes, collectProteins);
		}
	}

	/** Reject before even a default-model cold load can publish new process-wide geometry. */
	private static void requireLegacyCaller(){
		if(activeSession!=null){throw new IllegalStateException("Legacy caller cannot change settings during a batch session");}
	}

	/** Initializes resources before any worker starts; preserves the historical enabled callers. */
	private static void configureCaller(boolean boundaryNet){
		ProkObject.callCDS=ProkObject.calltRNA=ProkObject.call16S=true;
		ProkObject.call23S=ProkObject.call5S=true; ProkObject.call18S=false;
		ProkObject.loadSSUkmers=ProkObject.loadLSUkmers=ProkObject.load5Skmers=true;
		ProkObject.loadtRNAkmers=true;
		ProkObject.load16SSequence=ProkObject.load23SSequence=ProkObject.load5SSequence=true;
		CallGenes.loadTrnaResources(); ProkObject.loadLongKmers();
		ProkObject.loadConsensusSequenceFromFile(false, false);
		if(boundaryNet){loadDefaultTrnaBoundaryNet();}else{setBoundaryRefinement(false);}
	}

	/** Restores the same flags as the legacy single-bin adapter, including failed initialization. */
	private static final class CallerSettings implements AutoCloseable {
		final boolean oldCDS=ProkObject.callCDS, oldTRNA=ProkObject.calltRNA;
		final boolean old16S=ProkObject.call16S, old23S=ProkObject.call23S;
		final boolean old5S=ProkObject.call5S, old18S=ProkObject.call18S;
		final boolean oldLoadSSU=ProkObject.loadSSUkmers, oldLoadLSU=ProkObject.loadLSUkmers;
		final boolean oldLoad5S=ProkObject.load5Skmers, oldLoadTRNA=ProkObject.loadtRNAkmers;
		final boolean oldLoad16S=ProkObject.load16SSequence, oldLoad23S=ProkObject.load23SSequence;
		final boolean oldLoad5SSeq=ProkObject.load5SSequence, oldLoad18SSeq=ProkObject.load18SSequence;
		final boolean oldBoundary=boundaryRefinement();
		@Override public void close(){
			ProkObject.callCDS=oldCDS; ProkObject.calltRNA=oldTRNA;
			ProkObject.call16S=old16S; ProkObject.call23S=old23S;
			ProkObject.call5S=old5S; ProkObject.call18S=old18S;
			ProkObject.loadSSUkmers=oldLoadSSU; ProkObject.loadLSUkmers=oldLoadLSU;
			ProkObject.load5Skmers=oldLoad5S; ProkObject.loadtRNAkmers=oldLoadTRNA;
			ProkObject.load16SSequence=oldLoad16S; ProkObject.load23SSequence=oldLoad23S;
			ProkObject.load5SSequence=oldLoad5SSeq; ProkObject.load18SSequence=oldLoad18SSeq;
			setBoundaryRefinement(oldBoundary);
		}
	}

	/** Calls one bin using fixed session settings and a new bin-local caller/refinement model. */
	private static CallerData callConfiguredContigs(final ArrayList<Read> contigs, final GeneModel pgm,
			final int passes, final boolean collectProteins){
		assert(pgm!=null && passes>0) : "Caller geometry and refinement count must be validated before bin work";
		final GeneModel activePgm=CallGenes.makeMultipassModel(pgm, contigs, null, passes);
		final GeneCaller caller=CallGenes.makeGeneCaller(activePgm);
		final HashMap<String, CacheBuilder.Acc> accs=new HashMap<String, CacheBuilder.Acc>(contigs.size()*2);
		final ArrayList<ProteinSequence> proteins=(collectProteins ? new ArrayList<ProteinSequence>() : null);
		final HashMap<String, String> proteinToContig=(collectProteins ? new HashMap<String, String>() : null);
		for(final Read r : contigs){
			if(r==null || r.bases==null){continue;}
			final String name=canonicalName(r.id);
			if(accs.put(name, new CacheBuilder.Acc())!=null){
				throw new RuntimeException("Duplicate contig id: "+name);
			}
			final ArrayList<Orf> path=caller.callGenes(r, activePgm, true);
			final CacheBuilder.Acc acc=accs.get(name);
			if(path!=null){
				for(final Orf orf : path){addOrfChecked(acc, orf);}
				if(collectProteins){
					final ArrayList<Orf> cds=new ArrayList<Orf>();
					for(final Orf orf : path){if(orf.type==ProkObject.CDS){cds.add(orf);}}
					final ArrayList<Read> translated=CallGenes.translate(r, cds);
					if(translated!=null){for(final Read protein : translated){
						if(protein==null || protein.bases==null || protein.bases.length==0){continue;}
						final String proteinId=protein.id.replace('\t', '_');
						if(proteinToContig.put(proteinId, name)!=null){
							throw new RuntimeException("Duplicate translated protein id: "+proteinId);
						}
						proteins.add(new ProteinSequence(proteinId, protein.bases));
					}}
				}
			}
		}
		return new CallerData(accs, proteins, proteinToContig);
	}

	/**
	 * Opens one dedicated-JVM batch after taxonomy and all required PGMs have been loaded.
	 * Unrelated GeneCallAdapter/CallGenes invocations must not share this configured JVM.
	 */
	static synchronized CallerSession openSession(GeneModel[] models){
		if(activeSession!=null){throw new IllegalStateException("A caller batch session is already active");}
		if(models==null || models.length==0 || models[0]==null){throw new IllegalArgumentException("Preload the batch gene models first");}
		for(GeneModel model:models){
			if(!models[0].hasSameGeometry(model)){throw new IllegalArgumentException("Batch gene models have incompatible geometry");}
		}
		if(!models[0].matchesCurrentGeometry()){
			throw new IllegalStateException("Finish loading compatible gene models before opening the caller session");
		}
		final CallerSettings previous=new CallerSettings();
		try{
			configureCaller(true);
			activeSession=new CallerSession(models.clone(), previous);
			return activeSession;
		}catch(RuntimeException|Error failure){
			try{previous.close();}catch(RuntimeException|Error cleanup){failure.addSuppressed(cleanup);}
			throw failure;
		}
	}

	/** Owns setup lifetime; only entry/exit bookkeeping is serialized, never the per-bin work. */
	static final class CallerSession implements AutoCloseable {
		private CallerSession(GeneModel[] models_, CallerSettings previous_){models=models_; previous=previous_;}
		/** A worker uses only a preloaded parent and private searcher; refinement models stay local. */
		MagQCPreparedBin build(ArrayList<Read> contigs, String id, String status, String domain, String phylum,
				GeneModel model, int passes, ProteinSearcher.AssignmentBinding binding,
				ProteinSearcher searcher, ProteinSearcher.AssignPolicy policy, int lookahead){
			synchronized(FastaInCacheRowBuilder.class){
				if(activeSession!=this){throw new IllegalStateException("Caller session is closed");}
				boolean known=false;
				for(GeneModel candidate:models){if(candidate==model){known=true; break;}}
				if(!known){throw new IllegalArgumentException("Caller model was not preloaded for this batch");}
				active++; peakActive=Math.max(peakActive, active);
			}
			try{return buildPreparedBin(contigs, id, status, domain, phylum, model, passes,
				binding, searcher, policy, lookahead, true);}
			finally{synchronized(FastaInCacheRowBuilder.class){active--;}}
		}
		/** Maximum overlapping bin calls, useful for checking that a batch did not serialize its work. */
		int peakWorkers(){synchronized(FastaInCacheRowBuilder.class){return peakActive;}}
		/** Caller must join all workers before restoring the process-wide calling configuration. */
		@Override public void close(){
			synchronized(FastaInCacheRowBuilder.class){
				if(activeSession!=this){return;}
				if(active!=0){throw new IllegalStateException("Join caller workers before closing the batch session");}
				try{previous.close();}finally{activeSession=null;}
			}
		}
		private final GeneModel[] models;
		private final CallerSettings previous;
		private int active;
		private int peakActive;
	}

	/** Accessed only under the builder class monitor; legacy adapters cannot reconfigure an active batch. */
	private static CallerSession activeSession;

	/** Shared caller outputs; cache formatting and prepared aggregation consume the same accounting. */
	private static final class CallerData{
		CallerData(final HashMap<String,CacheBuilder.Acc> accs_,
				final ArrayList<ProteinSequence> proteins_,
				final HashMap<String,String> proteinToContig_){
			accs=accs_; proteins=proteins_; proteinToContig=proteinToContig_;
		}
		final HashMap<String,CacheBuilder.Acc> accs;
		final ArrayList<ProteinSequence> proteins;
		final HashMap<String,String> proteinToContig;
	}

	/** Checks the integer CDS coding accumulator before invoking the canonical Acc accounting. */
	private static void addOrfChecked(final CacheBuilder.Acc acc, final Orf orf){
		if(orf==null){return;}
		if(orf.type==ProkObject.CDS){
			long len=(long)orf.stop-(long)orf.start+1L;
			if(len<0){len=-len;}
			if(len<0 || len>Integer.MAX_VALUE || acc.cds==Integer.MAX_VALUE ||
				acc.coding>Integer.MAX_VALUE-len || acc.glenSum>Long.MAX_VALUE-len ||
				(len>0 && (len>Long.MAX_VALUE/len || acc.glenSq>Long.MAX_VALUE-len*len))){
				throw overflow("CDS accounting", "contig="+orf.scafName);
			}
		}else if(orf.type==ProkObject.tRNA && acc.trna==Integer.MAX_VALUE){
			throw overflow("trna_total", "contig="+orf.scafName);
		}else if(orf.type==ProkObject.r16S && acc.r16==Integer.MAX_VALUE){
			throw overflow("r16", "contig="+orf.scafName);
		}else if(orf.type==ProkObject.r23S && acc.r23==Integer.MAX_VALUE){
			throw overflow("r23", "contig="+orf.scafName);
		}else if(orf.type==ProkObject.r5S && acc.r5==Integer.MAX_VALUE){
			throw overflow("r5", "contig="+orf.scafName);
		}else if(orf.type==ProkObject.r18S && acc.rother==Integer.MAX_VALUE){
			throw overflow("rother", "contig="+orf.scafName);
		}
		acc.addOrf(orf);
	}

	/** Adds a nonnegative long with an explicit overflow guard. */
	private static long addPositive(final long old, final long add, final String field, final String id){
		if(add<0 || old<0 || old>Long.MAX_VALUE-add){throw overflow(field,id);} return old+add;
	}
	/** Adds a nonnegative int with an explicit overflow guard. */
	private static int addInt(final int old, final int add, final String field, final String id){
		if(add<0 || old<0 || old>Integer.MAX_VALUE-add){throw overflow(field,id);} return old+add;
	}
	/** Creates the uniform prepared-aggregation overflow exception. */
	private static RuntimeException overflow(final String field, final String id){return new ArithmeticException("Prepared "+field+" overflow for "+id);}

	/**
	 * Assigns each translated CDS to zero-or-one family via the canonical {@link
	 * ProteinSearcher#assignFamily} and records one family copy per ASSIGNED CDS at its roster rank
	 * ({@link ProteinSearcher.FamilyAssignment#familyIdx}). No fallback: a rejection records nothing.
	 * One {@link ProteinSearcher.ShortlistScratch} is allocated via {@link
	 * ProteinSearcher.AssignmentBinding#newScratch()} and reused across every CDS (this runs inside the
	 * {@code synchronized}, single-threaded {@link #buildRows}, matching assignFamily's own reuse
	 * contract). assignFamily's empty-Q guard returns {@code NO_VALID_KMER}; the generic rejection path
	 * below records nothing. This implements the accepted discard policy for CDSs with no valid 5-mer.
	 */
	private static void assignFamiliesViaBinding(final HashMap<String, CacheBuilder.Acc> accs,
			final ArrayList<ProteinSequence> proteins, final HashMap<String, String> proteinToContig,
			final ProteinSearcher.AssignmentBinding binding, final ProteinSearcher searcher,
			final ProteinSearcher.AssignPolicy policy, final int lookahead){
		assert(proteins!=null && proteinToContig!=null && searcher!=null)
			: "assignFamiliesViaBinding requires the binding path's non-null protein/contig maps + searcher "
			+"(buildRowsWithAssignments guarantees a non-null searcher, so buildRows builds both maps at 148-149).";
		final ProteinSearcher.ShortlistScratch scratch=binding.newScratch();
		for(final ProteinSequence protein : proteins){
			final ProteinSearcher.FamilyAssignment result=
				searcher.assignFamily(binding,protein,scratch,policy,lookahead);
			if(!result.isAssigned()){continue;}//rejection -> no family copy, NO fallback
			final String contig=proteinToContig.get(protein.id);
			if(contig==null){throw new RuntimeException("Missing contig owner for translated protein "+protein.id);}
			final CacheBuilder.Acc acc=accs.get(contig);
			if(acc==null){throw new RuntimeException("Missing accumulator for translated protein contig "+contig);}
			final int old=acc.fam.get(result.familyIdx);
			if(old==Integer.MAX_VALUE || acc.mapped==Integer.MAX_VALUE){
				throw overflow("family assignment", contig);
			}
			acc.fam.put(result.familyIdx, old<0 ? 1 : old+1);
			acc.mapped++;
		}
	}

	/**
	 * Mirrors CallGenes.parse()'s bare-invocation boundary-net default.  The
	 * BBTools method is package-private in {@code prok.TrnaCaller}; reflection is
	 * confined to this compatibility bridge so the public BBTools tree remains
	 * untouched.  An already configured boundary state is preserved.
	 */
	private static void loadDefaultTrnaBoundaryNet(){
		try{
			final Class<?> clazz=Class.forName("prok.TrnaCaller");
			final Field refine=boundaryField(clazz);
			if(refine.getBoolean(null)){return;}
			final String net=Data.findPath("?trna_boundary_net.bbnet");
			final String start=Data.findPath("?trna_boundary_start_table.tsv");
			final String stop=Data.findPath("?trna_boundary_stop_table.tsv");
			if(net==null || start==null || stop==null){
				throw new RuntimeException("CallGenes boundary resources unresolved: net="+net+
					" start="+start+" stop="+stop);
			}
			final Method load=clazz.getDeclaredMethod("loadBoundaryNet", String.class, String.class,
				String.class, String.class);
			load.setAccessible(true);
			load.invoke(null, net, net, start, stop);
		}catch(final RuntimeException e){throw e;
		}catch(final Exception e){throw new RuntimeException("Could not mirror CallGenes boundary-net default", e);}
	}

	/** Reads the package-private BBTools boundary state through the same bridge. */
	private static boolean boundaryRefinement(){
		try{return boundaryField(Class.forName("prok.TrnaCaller")).getBoolean(null);}
		catch(final Exception e){throw new RuntimeException("Could not inspect CallGenes boundary state", e);}
	}

	/** Sets the package-private BBTools boundary state for the focused regression. */
	private static void setBoundaryRefinement(final boolean value){
		try{boundaryField(Class.forName("prok.TrnaCaller")).setBoolean(null, value);}
		catch(final Exception e){throw new RuntimeException("Could not restore CallGenes boundary state", e);}
	}

	/** Resolves the package-private state field without changing the BBTools source tree. */
	private static Field boundaryField(final Class<?> clazz) throws Exception{
		final Field refine=clazz.getDeclaredField("REFINE_BOUNDARIES");
		refine.setAccessible(true);
		return refine;
	}

	/** CacheBuilder's canonical key after the training path's trd=t normalization. */
	private static String canonicalName(final String id){
		if(id==null || id.length()==0){throw new RuntimeException("Empty contig id.");}
		int end=id.length();
		for(int i=0; i<end; i++){
			if(Character.isWhitespace(id.charAt(i))){end=i; break;}
		}
		final StringBuilder sb=new StringBuilder(end);
		for(int i=0; i<end; i++){sb.append(id.charAt(i)==' ' ? '_' : id.charAt(i));}
		return sb.toString();
	}
}
