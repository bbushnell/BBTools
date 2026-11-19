package prok;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;
import java.util.regex.Pattern;

import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaEvidenceCatalog.Membership;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Assembly-only holdout audit before construction. Exact short-gene identity is
 * deliberately absent. Canonical file aliases and explicit NCBI assembly-record
 * identities establish exclusion groups. Equal genome dimensions only nominate
 * possible duplicate assemblies for subsequent sequence comparison; they never
 * establish an alias. No sequence payload, extraction or training is performed.
 * @author Ganyu
 */
public final class NcrnaConstructionAssemblies{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaConstructionAssemblies.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String original=null, ec=null, em=null, dc=null, dm=null;
			int oldCount=-1, evalCount=-1, newCount=-1;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("genomes")){original=b;}
				else if(a.equals("evalcensus")){ec=b;}
				else if(a.equals("evalmembers")){em=b;}
				else if(a.equals("newcensus")){dc=b;}
				else if(a.equals("newmembers")){dm=b;}
				else if(a.equals("expectedgenomes")){oldCount=Parse.parseIntKMG(b);}
				else if(a.equals("expectedeval")){evalCount=Parse.parseIntKMG(b);}
				else if(a.equals("expectednew")){newCount=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown assembly-audit parameter: "+arg);}
			}
			boolean assertions=false;
			assert(assertions=true);
			require(assertions && original!=null && ec!=null && em!=null && dc!=null && dm!=null && parser.out1!=null,
				"Supply -ea genomes= evalcensus= evalmembers= newcensus= newmembers= and out=fresh_prefix");
			require(parser.in1==null && parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
				&& parser.maxReads<0 && parser.samplerate==1, "Assembly audit must retain the complete declared cohorts");
			process(original, NcrnaEvidenceCatalog.load(ec, em), NcrnaEvidenceCatalog.load(dc, dm),
				oldCount, evalCount, newCount, parser.out1);
			pp.outstream.println("NCRNA_CONSTRUCTION_ASSEMBLY_AUDIT_COMPLETE_NOT_FROZEN_SPLIT");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public static void process(final String original, final NcrnaEvidenceCatalog evaluation,
			final NcrnaEvidenceCatalog discovery, final int oldCount, final int evalCount, final int newCount, final String output){
		require(original!=null && evaluation!=null && discovery!=null && output!=null && oldCount>0
			&& evalCount>0 && newCount>0 && evaluation.size()==evalCount && discovery.size()==newCount,
			"Assembly audit needs complete original/evaluation/discovery population counts");
		RrnaResourceIO.fresh(output, ".inputs.tsv.gz", ".candidate_pairs.tsv.gz", ".summary.tsv");
		final NcrnaConstructionAssemblies audit=new NcrnaConstructionAssemblies();
		audit.original(original, oldCount);
		audit.catalog(evaluation, true);
		audit.catalog(discovery, false);
		audit.write(output);
	}

	private NcrnaConstructionAssemblies(){}

	/*--------------------------------------------------------------*/
	/*----------------        Assembly Binding      ----------------*/
	/*--------------------------------------------------------------*/

	private void original(final String path, final int expected){
		int count=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionBindings.GENOMES_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==11 && p.parseInt(0)==count++, "Original assembly roster must preserve its complete ordinal order");
				p.setBounds(7);
				final String fasta=new String(row, p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
				final String assembly=p.parseString(2), real=canonical(fasta);
				final long loci=p.parseLong(8), bases=p.parseLong(9), records=p.parseLong(10);
				// NcrnaAnnotatedCorpus.extractGenome returns before reading FASTA
				// when loci is empty; its paired zeros mean unmeasured dimensions.
				require(accession(assembly)!=null && (loci>0 && bases>0 && records>0 || loci==0 && bases==0 && records==0),
					"Original dimensions must be measured for target-bearing genomes or paired unmeasured zeros for zero loci");
				if(bases==0){unmeasuredOriginal++;}
				final Original old=originals.get(real);
				if(old==null){originals.put(real, new Original(assembly, bases, records));}
				else{
					require(old.assembly.equals(assembly) && (old.bases==0 || bases==0 || old.bases==bases && old.records==records),
						"One canonical original file cannot carry conflicting measured assembly metadata");
					if(old.bases==0 && bases>0){originals.put(real, new Original(assembly, bases, records));}
				}
			}
		}
		require(count==expected, "Original genome roster count differs from its declared population");
	}

	private void catalog(final NcrnaEvidenceCatalog catalog, final boolean evaluation){
		assert(catalog.size()>0) : "Both cohort denominators were checked before assembly binding";
		for(int i=0; i<catalog.size(); i++){
			final Genome genome=catalog.get(i);
			final String real=canonical(genome.fasta);
			require(new File(real).length()==genome.bytes, "A catalog source file changed size since its saved census");
			final Node node=new Node(nodes.size(), evaluation, genome, real);
			nodes.add(node);
			final Original old=originals.get(real);
			if(old!=null){
				require(old.bases==0 || old.bases==genome.bases && old.records==genome.records,
					"Canonical original/catalog aliases must preserve dimensions actually observed by the original extractor");
				identity(node, old.assembly);
			}
			for(int m=0; m<genome.memberships(); m++){
				final Membership member=genome.membership(m);
				if(member.assembly.equals("NA")){continue;}
				// Membership records are source claims, not proof that a supplied
				// FASTA is the census file. A foreign/unavailable claim cannot assign
				// this input to an evaluation assembly, even if dimensions agree.
				if(!new File(member.suppliedFasta).isFile() || !real.equals(canonical(member.suppliedFasta))){
					node.unresolvedClaim=true;
				}else{identity(node, member.assembly);}
			}
			final Node sameFile=files.get(real);
			if(sameFile==null){files.put(real, node);}
			else{
				require(sameFile.genome.bases==genome.bases && sameFile.genome.records==genome.records,
					"The same canonical FASTA cannot have different observed nucleotide/record counts");
				union(node, sameFile);
			}
			if(node.assembly!=null){
				final Node sameAssembly=assemblies.get(node.assembly);
				if(sameAssembly==null){assemblies.put(node.assembly, node);}
				else{union(node, sameAssembly);}
			}
			final String dimensions=genome.bases+"/"+genome.records;
			Dimensions group=dimensionIndex.get(dimensions);
			if(group==null){
				group=new Dimensions();
				dimensionIndex.put(dimensions, group);
				dimensionOrder.add(group);
			}
			(evaluation ? group.evaluation : group.discovery).add(node);
		}
	}

	private static void identity(final Node node, final String claim){
		final String assembly=accession(claim);
		if(assembly==null){
			node.unresolvedClaim=true;
			return;
		}
		require(node.assembly==null || node.assembly.equals(assembly),
			"Conflicting assembly-record claims require explicit resolution, not an arbitrary chosen identity");
		node.assembly=assembly;
	}

	/** Versions belong to the same declared assembly record for holdout exclusion.
	 * GCA and GCF accessions are not equated without an explicit binding. */
	private static String accession(final String claim){
		if(claim==null || !ACCESSION.matcher(claim).matches()){return null;}
		final int dot=claim.indexOf('.');
		return dot<0 ? claim : claim.substring(0, dot);
	}

	private String canonical(final String path){
		require(path!=null && path.indexOf('\ufffd')<0, "Invalid UTF-8 must not silently alter assembly path identity");
		String result=canonicalPaths.get(path);
		if(result!=null){return result;}
		require(path!=null && new File(path).isAbsolute() && new File(path).isFile(),
			"Assembly aliases require existing absolute FASTA files, not guessed or unavailable paths");
		try{result=new File(path).toPath().toRealPath().toString();}
		catch(IOException error){throw new IllegalArgumentException("Could not resolve assembly source: "+path, error);}
		canonicalPaths.put(path, result);
		return result;
	}

	private Node root(final Node node){
		Node result=node;
		while(result.parent!=result.index){result=nodes.get(result.parent);}
		Node current=node;
		while(current.parent!=result.index){
			final Node next=nodes.get(current.parent);
			current.parent=result.index;
			current=next;
		}
		return result;
	}
	private void union(final Node first, final Node second){
		final Node a=root(first), b=root(second);
		if(a==b){return;}
		final Node keep=a.index<b.index ? a : b, drop=keep==a ? b : a;
		drop.parent=keep.index;
		keep.evaluationGroup|=drop.evaluationGroup;
		keep.unresolvedClaim|=drop.unresolvedClaim;
		if(keep.assembly==null){keep.assembly=drop.assembly;}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Explicit Audit        ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String prefix){
		long pairs=0, excluded=0, unresolved=0, candidates=0, train=0, eval=0;
		try(Output inputs=new Output(prefix+".inputs.tsv.gz"); Output comparisons=new Output(prefix+".candidate_pairs.tsv.gz");
				Output summary=new Output(prefix+".summary.tsv")){
			comparisons.writer.println("evaluation_input\tdiscovery_input\tevaluation_fasta\tdiscovery_fasta\tbases\tcontigs\tbasis");
			for(final Dimensions group:dimensionOrder){
				for(final Node discovery:group.discovery){
					if(root(discovery).evaluationGroup){continue;}
					for(final Node evaluation:group.evaluation){
						root(discovery).sequenceCandidate=true;
						pairs++;
						require(pairs<=1000000, "More than one million candidate assembly pairs requires a revised bounded comparison design");
						comparisons.writer.print(buffer.clear().append(evaluation.genome.id).tab().append(discovery.genome.id).tab()
							.append(evaluation.fastaBytes).tab().append(discovery.fastaBytes).tab().append(discovery.genome.bases)
							.tab().append(discovery.genome.records).append("\tDIMENSION_COMPATIBILITY_NOT_SEQUENCE_IDENTITY\n"));
					}
				}
			}
			inputs.writer.println("cohort\tinput_id\tfasta\tcanonical_fasta\tassembly_record\tgroup_id\tevaluation_group\tsequence_comparison_pending\tunresolved_claim\tdisposition");
			for(final Node node:nodes){
				final Node group=root(node);
				final String disposition;
				if(node.evaluation){
					eval++;
					disposition="EVALUATION";
				}else if(group.evaluationGroup){
					excluded++;
					disposition="EXCLUDE_KNOWN_EVALUATION_ASSEMBLY";
				}else if(group.sequenceCandidate){
					candidates++;
					disposition="ASSEMBLY_SEQUENCE_CHECK_PENDING";
				}else if(group.assembly==null || group.unresolvedClaim){
					unresolved++;
					disposition="ASSEMBLY_IDENTITY_UNRESOLVED";
				}else{
					train++;
					disposition="TRAIN_CANDIDATE_ASSEMBLY_IDENTITY";
				}
				inputs.writer.print(buffer.clear().append(node.evaluation ? "evaluation" : "discovery").tab().append(node.genome.id)
					.tab().append(node.fastaBytes).tab().append(node.canonicalBytes).tab().append(node.assembly==null ? "NA" : node.assembly)
					.tab().append(group.index).tab().append(group.evaluationGroup).tab().append(group.sequenceCandidate)
					.tab().append(group.unresolvedClaim).tab().append(disposition).nl());
			}
			assert(eval+excluded+unresolved+candidates+train==nodes.size()) : "Disjoint assembly dispositions must conserve every input in both cohorts";
			summary.writer.println("evaluation_inputs\tknown_excluded_discovery\tsequence_check_discovery\tunresolved_discovery\ttraining_candidate_discovery\tcandidate_pairs\tunmeasured_original_genomes\tstatus");
			summary.writer.print(buffer.clear().append(eval).tab().append(excluded).tab().append(candidates).tab().append(unresolved)
				.tab().append(train).tab().append(pairs).tab().append(unmeasuredOriginal).append("\tASSEMBLY_AUDIT_NOT_FROZEN_MEMBER_SPLIT\n"));
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            State             ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Original{
		Original(final String a, final long b, final long r){
			assembly=a;
			bases=b;
			records=r;
		}
		final String assembly;
		final long bases, records;
	}
	private static final class Node{
		Node(final int index_, final boolean evaluation_, final Genome genome_, final String canonical_){
			index=index_;
			parent=index;
			evaluation=evaluation_;
			evaluationGroup=evaluation;
			genome=genome_;
			canonical=canonical_;
			fastaBytes=genome.fasta.getBytes(StandardCharsets.UTF_8);
			canonicalBytes=canonical.getBytes(StandardCharsets.UTF_8);
		}
		final int index;
		final boolean evaluation;
		final Genome genome;
		final String canonical;
		final byte[] fastaBytes, canonicalBytes;
		int parent;
		String assembly;
		boolean evaluationGroup, unresolvedClaim, sequenceCandidate;
	}
	private static final class Dimensions{
		final ArrayList<Node> evaluation=new ArrayList<Node>(), discovery=new ArrayList<Node>();
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private final ObjectMap<String, String> canonicalPaths=new ObjectMap<String, String>(String.class, String.class);
	private final ObjectMap<String, Original> originals=new ObjectMap<String, Original>(String.class, Original.class);
	private final ObjectMap<String, Node> files=new ObjectMap<String, Node>(String.class, Node.class);
	private final ObjectMap<String, Node> assemblies=new ObjectMap<String, Node>(String.class, Node.class);
	private final ObjectMap<String, Dimensions> dimensionIndex=new ObjectMap<String, Dimensions>(String.class, Dimensions.class);
	private final ArrayList<Dimensions> dimensionOrder=new ArrayList<Dimensions>();
	private final ArrayList<Node> nodes=new ArrayList<Node>();
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder();
	private long unmeasuredOriginal;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private static final Pattern ACCESSION=Pattern.compile("GC[AF]_[0-9]+(?:\\.[0-9]+)?");
}
