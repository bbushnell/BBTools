package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import map.IntHashMap;
import map.LongIntMap;
import map.ObjectIntMap;
import parse.LineParser1;
import prok.GradeNcrnaReference.Binding;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/** Describes saved source-to-event-to-parent correspondence, never re-adjudicates.
 * Exact selected source intervals collapse within input/contig/family/strand.
 * Only valid same-family links enter the scored bipartite graph. Unassigned
 * sources, masks, protected priors and missing paths remain explicit. Connected
 * components, not arbitrary best-overlap matches, define splits and merges.
 * Weak/bridge witness associations remain diagnostic paths outside components.
 * One thread owns primitive indexes and reusable byte probes.
 * @author Ganyu
 */
final class NcrnaSourceComparison{

	static void validatePaths(final String directory, final String reference, final String prefix){
		final ArrayList<String> paths=new ArrayList<String>();
		paths.add(reference);
		for(final String name:INPUTS){paths.add(new File(directory, name).getPath());}
		require(Tools.testInputFiles(false, true, paths.toArray(new String[paths.size()])), "Saved comparison inputs must be readable");
		for(final String suffix:SUFFIXES){
			require(Tools.testOutputFiles(false, false, false, prefix+suffix), "Comparison output must be fresh: "+prefix+suffix);
			paths.add(prefix+suffix);
		}
		require(Tools.testForDuplicateFiles(true, paths.toArray(new String[paths.size()])), "Comparison input/output paths must be distinct");
	}

	static void process(final String directory, final String reference, final String prefix){
		final NcrnaSourceComparison report=new NcrnaSourceComparison(directory, prefix);
		report.readParents(reference);
		report.readEvents();
		report.readWitnesses();
		report.readSources();
		report.readLinks();
		report.checkCounts();
		report.components();
		report.writeResults();
	}

	private NcrnaSourceComparison(final String directory_, final String prefix_){
		directory=directory_;
		prefix=prefix_;
	}

	/*----------------       Saved populations       ----------------*/

	private void readParents(final String reference){
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		for(final Binding binding:GradeNcrnaReference.loadManifest(reference)){
			final int input=symbol(binding.id);
			require(!inputCounts.contains(input), "Repeated reference input identity");
			inputCounts.put(input, inputs.size());
			final Population population=new Population(input);
			inputs.add(population);
			for(final Row row:codec.read(binding.gff, binding.id)){
				final Parent parent=new Parent(input, symbol(row.locus.seqid), symbol(row.locus.family),
					row.locus.start, row.locus.stop, row.locus.strand, symbol(row.id), row.mask);
				key.clear().append(input).tab().append(parent.name);
				require(parentKeys.get(probe(key))<0, "Repeated input/parent identity");
				parentKeys.put(probe(key).copy(), parents.size());
				parents.add(parent);
				population.counts[row.mask ? 1 : 0]++;
				family(parent.family).observed[row.mask ? 9 : 8]++;
			}
		}
		require(!inputs.isEmpty(), "Comparison requires a nonempty reference manifest, including empty-reference inputs");
	}

	private void readEvents(){
		try(Input in=open("event_provenance.tsv.gz", EVENT_HEADER)){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parse(line, 30);
				final int id=positive(0);
				final Event event=new Event(input(1), symbol(2), symbol(7), parser.parseInt(3), parser.parseInt(4), strand(5),
					symbol(6), symbol(26), parser.parseLong(27), state(EVENT_STATES, 22));
				require(event.line>0 && eventIds.get(id)<0, "Native event IDs must be unique with positive native line");
				eventIds.put(id, events.size());
				events.add(event);
				final Family family=family(event.family);
				family.observed[0]++;
				family.observed[event.tier+1]++;
			}
		}
	}

	private void readWitnesses(){
		try(Input in=open("annotation_witnesses.tsv.gz", WITNESS_HEADER)){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parse(line, 11);
				final int input=input(1);
				key.clear().append(input).tab().append(symbol(0));
				final int index=parentKeys.get(probe(key));
				require(index>=0, "Annotation witness has no final parent in the declared reference");
				final Parent parent=parents.get(index);
				if(parser.termEquals("PROTECTED_PRIOR", 2)){
					require(!parent.prior && !parent.mask && parser.termEquals(".", 3)
						&& parser.termEquals(".", 4) && parser.termEquals(".", 5) && parser.parseLong(9)>0
						&& parser.parseInt(7)==parent.start && parser.parseInt(8)==parent.stop
						&& parser.termEquals("ORIGINAL_ACCEPTED_INTERVAL_AND_AUTHORITY", 10), "Malformed or duplicate protected-prior witness");
					parent.prior=true;
					inputs.get(inputCounts.get(input)).counts[2]++;
					family(parent.family).observed[10]++;
				}else{
					require(parser.termEquals("CM_HIT", 2), "Unknown annotation evidence kind");
					final Event event=event(positive(3));
					require(event.input==input && event.target==parent.target && event.family==parent.family
						&& event.strand==parent.strand && event.table==symbol(4) && event.line==parser.parseLong(5)
						&& event.model==symbol(6) && event.start==parser.parseInt(7) && event.stop==parser.parseInt(8),
						"Annotation witness must bind the original event identity, geometry and native source");
					final int association=state(ASSOCIATIONS, 10);
					require(!contains(event.parents, index), "Repeated parent/event witness would duplicate a correspondence path");
					if(event.parents==null){
						event.parents=new IntList(1);
						event.associations=new IntList(1);
					}
					event.parents.add(index);
					event.associations.add(association);
					parent.witnesses++;
					if(association==0){parent.priorAssociations++;}
				}
			}
		}
		for(final Parent parent:parents){
			require(parent.prior || parent.witnesses>0, "Every final annotation needs a retained prior or CM witness");
			require(parent.priorAssociations==0 || parent.prior, "Prior-associated CM witness needs a protected-prior witness");
		}
	}

	/** Copies selected claims to the alias ledger. Positive indexes name selected
	 * rows; -(family+2) retains a context row's family, while -1 means absent. */
	private void readSources(){
		try(Input in=open("source_provenance.tsv.gz", SOURCE_HEADER); Output out=new Output(prefix+".source_rows.tsv.gz")){
			out.writer.println("physical_source_id\t"+SOURCE_HEADER);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parse(line, 13);
				final int input=input(0);
				final int file=symbol(1);
				final long lineNumber=parser.parseLong(2);
				final int family=symbol(8);
				final int disposition=state(SOURCE_STATES, 10);
				final boolean selected=(parser.parseInt(9)&3)!=0;
				require(lineNumber>0 && (disposition!=0)==selected, "Selected-source bits must agree with saved disposition");
				family(family).sourceCounts[disposition]++;
				key.clear().append(input).tab().append(file);
				int fileIndex=fileKeys.get(probe(key));
				if(fileIndex<0){
					fileIndex=files.size();
					fileKeys.put(probe(key).copy(), fileIndex);
					files.add(new LongIntMap());
				}
				final LongIntMap rows=files.get(fileIndex);
				require(!rows.contains(lineNumber), "Repeated original source input/file/line key");
				int index=-1;
				if(selected){
					final int target=symbol(3);
					final int start=parser.parseInt(5);
					final int stop=parser.parseInt(6);
					final char strand=strand(7);
					key.clear().append(input).tab().append(target).tab().append(family).tab().append(start).tab().append(stop).tab().append(strand);
					index=sourceKeys.get(probe(key));
					if(index<0){
						index=sources.size();
						sourceKeys.put(probe(key).copy(), index);
						sources.add(new Source(input, target, family, start, stop, strand));
					}
					final Source source=sources.get(index);
					source.aliases++;
					source.dispositions|=1<<disposition;
					out.writer.print(buffer.clear().append(index+1).tab().append(line).nl());
				}
				rows.put(lineNumber, selected ? index+1 : -family-2);
			}
		}
		sourceKeys.clear();
	}

	/** Retains every qualifying path in a ledger but counts a physical graph edge
	 * once even if duplicated source claims or multiple models support it. */
	private void readLinks(){
		try(Input in=open("source_event_links.tsv.gz", LINK_HEADER); Output out=new Output(prefix+".paths.tsv.gz")){
			out.writer.println("physical_source_id\tannotation_id\tsource_input\tsource_gff\tsource_line\tevent_id\tparent_mask\tfamily_agreement\tassociation");
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parse(line, 7);
				final int input=input(0);
				key.clear().append(input).tab().append(symbol(1));
				final int file=fileKeys.get(probe(key));
				require(file>=0, "Source link names an absent input/file");
				final int stored=files.get(file).get(parser.parseLong(2));
				require(stored!=-1 && stored!=0, "Source link names an absent original source line");
				final int eventId=positive(3);
				final Event event=event(eventId);
				require(event.input==input && parser.termEquals(EVENT_STATES[event.tier], 6), "Source link and event identity/tier disagree");
				final boolean valid=bool(4);
				final int agreement=state(AGREEMENTS, 5);
				if(stored<-1){
					contextLinks++;
					family(-stored-2).contextLinks++;
					continue;
				}
				final int index=stored-1;
				final Source source=sources.get(index);
				require(source.target==event.target, "Saved source link crosses contigs");
				final int expected=isUnassigned(source) ? 2 : source.family==event.family ? 0 : 1;
				require(agreement==expected, "Source link family agreement differs from its declared source and event families");
				source.eventLinks++;
				if(!valid){continue;}
				require(source.start>0 && source.stop>=source.start && source.strand==event.strand
					&& (source.strand=='+' || source.strand=='-'), "Valid source link requires assessable, strand-consistent geometry");
				if(event.parents==null){source.unparentedEvents++; continue;}
				for(int j=0; j<event.parents.size; j++){
					final int parentIndex=event.parents.get(j);
					final Parent parent=parents.get(parentIndex);
					buffer.clear().append(index+1).tab().append(name(parent.name)).tab();
					appendTerm(buffer, 0).tab();
					appendTerm(buffer, 1).tab().append(parser.parseLong(2)).tab().append(eventId).tab().append(parent.mask).tab()
						.append(AGREEMENTS[agreement]).tab().append(ASSOCIATIONS[event.associations.get(j)]).nl();
					out.writer.print(buffer);
					final int association=event.associations.get(j);
					if(association==1 || association==2){
						family(source.family).diagnosticPaths[association-1]++;
						if(!parent.mask){parent.diagnosticSources=add(parent.diagnosticSources, index);}
						if(agreement!=1 && !parent.mask){
							source.diagnosticScored=add(source.diagnosticScored, parentIndex);
						}
						continue;
					}
					require(association==0 || association==3, "Only explicitly declared prior or qualifying-group associations can enter correspondence components");
					if(agreement!=1 && !parent.mask){source.directScored=true;}
					if(agreement==1){
						parent.otherFamilies=add(parent.otherFamilies, index);
						source.otherParents=add(source.otherParents, parentIndex);
					}else if(agreement==2){
						parent.unassigned=add(parent.unassigned, index);
						source.unassignedParents=add(source.unassignedParents, parentIndex);
					}else if(parent.mask){
						source.masks=add(source.masks, parentIndex);
						parent.sources=add(parent.sources, index);
					}else{
						addDirect(source, parentIndex, association==0 ? PRIOR_KIND : GROUP_KIND);
						parent.sources=add(parent.sources, index);
					}
				}
			}
		}
		files.clear();
		fileKeys.clear();
	}

	/*----------------       Conservation gates       ----------------*/

	private void checkCounts(){
		try(Input in=open("family_counts.tsv", FAMILY_HEADER)){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parse(line, 12);
				final Family family=family(symbol(0));
				require(!family.checked, "Repeated saved family accounting");
				family.checked=true;
				for(int i=1; i<12; i++){
					require(parser.parseLong(i)>=0 && parser.parseLong(i)==family.observed[i-1], "Event/parent counts differ from saved family accounting: "+name(family.name));
				}
			}
		}
		try(Input in=open("source_family_counts.tsv", SOURCE_COUNT_HEADER)){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parse(line, 13);
				final Family family=family(symbol(0));
				require(!family.sourceChecked, "Repeated saved source-family accounting");
				family.sourceChecked=true;
				long sum=0;
				for(int i=0; i<SOURCE_STATES.length; i++){
					require(parser.parseLong(i+2)==family.sourceCounts[i], "Source rows differ from saved disposition totals: "+name(family.name));
					sum+=family.sourceCounts[i];
				}
				require(parser.parseLong(1)==sum, "Source disposition totals do not conserve their saved total");
			}
		}
		try(Input in=open("summary.tsv", INPUT_HEADER)){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parse(line, 5);
				final Population pop=inputs.get(inputCounts.get(input(0)));
				require(!pop.checked, "Repeated saved input summary");
				pop.checked=true;
				for(int i=0; i<3; i++){require(parser.parseLong(i+1)==pop.counts[i], "Per-input scored/masked/prior totals disagree");}
			}
		}
		for(final Population pop:inputs){require(pop.checked, "Reference manifest contains an input absent from saved summary");}
		for(final Family family:families){
			require(family.checked || sum(family.observed)==0, "Observed events/parents have no saved family accounting");
			require(family.sourceChecked || sum(family.sourceCounts)==0, "Observed source rows have no saved source accounting");
		}
	}

	/*----------------       Physical components       ----------------*/

	private void components(){
		final IntList sourceQueue=new IntList();
		final IntList parentQueue=new IntList();
		try(Output out=new Output(prefix+".components.tsv.gz")){
			out.writer.println("component_id\tfamily\tsource_loci\tscored_parents\tshape");
			for(int i=0; i<sources.size(); i++){
				final Source first=sources.get(i);
				if(first.component!=0 || size(first.parents)==0){continue;}
				components++;
				sourceQueue.clear();
				parentQueue.clear();
				sourceQueue.add(i);
				first.component=components;
				for(int next=0; next<sourceQueue.size; next++){
					final Source source=sources.get(sourceQueue.get(next));
					assert(source.family==first.family) : "Only same-family saved links belong in a scored correspondence component";
					for(int j=0; j<source.parents.size; j++){
						final int p=source.parents.get(j);
						final Parent parent=parents.get(p);
						if(parent.component!=0){continue;}
						parent.component=components;
						parentQueue.add(p);
						for(int k=0; k<parent.sources.size; k++){
							final int s=parent.sources.get(k);
							if(sources.get(s).component==0){
								sources.get(s).component=components;
								sourceQueue.add(s);
							}
						}
					}
				}
				final int shape=sourceQueue.size==1 ? (parentQueue.size==1 ? PAIR : SPLIT) : parentQueue.size==1 ? MERGE : COMPLEX;
				for(int j=0; j<sourceQueue.size; j++){sources.get(sourceQueue.get(j)).shape=shape;}
				for(int j=0; j<parentQueue.size; j++){parents.get(parentQueue.get(j)).shape=shape;}
				out.writer.print(buffer.clear().append(components).tab().append(name(first.family)).tab()
					.append(sourceQueue.size).tab().append(parentQueue.size).tab().append(SHAPES[shape]).nl());
			}
		}
	}

	private void writeResults(){
		try(Output out=new Output(prefix+".sources.tsv.gz")){
			out.writer.println("physical_source_id\tinput_id\ttarget\tfamily\tstart\tstop\tstrand\tclaim_rows\tcomponent_id\tstatus\tscored_parents\tmasked_parents\tunassigned_parent_links\tother_family_parents\tweak_or_bridge_scored_parents\tevent_links\tunparented_events\tpair_association\tshift_5prime\tshift_3prime\tsource_length_minus_parent_length");
			for(int i=0; i<sources.size(); i++){
				final Source source=sources.get(i);
				final Family family=family(source.family);
				final String status=sourceStatus(source);
				family.sources++;
				family.statuses[state(SOURCE_RESULTS, status)]++;
				buffer.clear().append(i+1).tab();
				geometry(buffer, source).tab().append(source.aliases).tab().append(source.component).tab().append(status)
					.tab().append(size(source.parents)).tab().append(size(source.masks)).tab().append(size(source.unassignedParents))
					.tab().append(size(source.otherParents)).tab().append(size(source.diagnosticScored)).tab()
					.append(source.eventLinks).tab().append(source.unparentedEvents);
				if(source.shape==PAIR){
					final Parent parent=parents.get(source.parents.get(0));
					final int kind=source.parentKinds.get(0);
					final int shift5=source.strand=='+' ? parent.start-source.start : source.stop-parent.stop;
					final int shift3=source.strand=='+' ? parent.stop-source.stop : source.start-parent.start;
					final long length=(long)source.stop-source.start-((long)parent.stop-parent.start);
					buffer.tab().append(PAIR_KINDS[kind]).tab().append(shift5).tab().append(shift3).tab().append(length);
					if(kind==GROUP_KIND && (shift5!=0 || shift3!=0)){
						family.shift5.add(shift5);
						family.shift3.add(shift3);
						family.lengths[length<0 ? 0 : length==0 ? 1 : 2]++;
					}
				}else{buffer.append("\tNA\tNA\tNA\tNA");}
				out.writer.print(buffer.nl());
			}
		}
		try(Output out=new Output(prefix+".parents.tsv.gz")){
			out.writer.println("annotation_id\tinput_id\ttarget\tfamily\tstart\tstop\tstrand\tcomponent_id\tstatus\tprotected_prior\tcm_witnesses\tlinked_source_loci\tunassigned_source_loci\tother_family_source_loci\tweak_or_bridge_source_loci");
			for(final Parent parent:parents){
				final String status=parentStatus(parent);
				family(parent.family).parentStatuses[state(PARENT_RESULTS, status)]++;
				buffer.clear().append(name(parent.name)).tab();
				geometry(buffer, parent).tab().append(parent.component).tab().append(status).tab().append(parent.prior)
					.tab().append(parent.witnesses).tab().append(size(parent.sources)).tab().append(size(parent.unassigned))
					.tab().append(size(parent.otherFamilies)).tab().append(size(parent.diagnosticSources));
				out.writer.print(buffer.nl());
			}
		}
		writeShifts();
		// Publish completion accounting only after all detail writers close cleanly.
		try(Output out=new Output(prefix+".summary.tsv")){
			buffer.clear().append("family\tinputs\tselected_source_claims\tphysical_sources\texact_duplicate_claims\tcontext_source_rows\tcontext_links");
			for(final String status:SOURCE_RESULTS){buffer.tab().append(status);}
			for(final String status:PARENT_RESULTS){buffer.tab().append(status);}
			out.writer.print(buffer.append("\tsource_shorter_changed_pairs\tequal_length_changed_pairs\tsource_longer_changed_pairs\tbridge_paths\tweak_overlap_paths\tqualification\n"));
			long contextTotal=0;
			for(final Family family:families){
				final long selected=sum(family.sourceCounts)-family.sourceCounts[0];
				require(sum(family.statuses)==family.sources && selected>=family.sources, "Source classification must conserve collapsed loci");
				require(sum(family.parentStatuses)==family.observed[8]+family.observed[9], "Parent classification must conserve scored plus masked annotations");
				buffer.clear().append(name(family.name)).tab().append(inputs.size()).tab().append(selected).tab().append(family.sources)
					.tab().append(selected-family.sources).tab().append(family.sourceCounts[0]).tab().append(family.contextLinks);
				contextTotal+=family.contextLinks;
				for(final long n:family.statuses){buffer.tab().append(n);}
				for(final long n:family.parentStatuses){buffer.tab().append(n);}
				for(final long n:family.lengths){buffer.tab().append(n);}
				for(final long n:family.diagnosticPaths){buffer.tab().append(n);}
				out.writer.print(buffer.append("\tSAVED_LINK_CORRESPONDENCE_NOT_NEW_ADJUDICATION_OR_BIOLOGICAL_ABSENCE\n"));
			}
			require(contextTotal==contextLinks, "Declared source-family context links must conserve the total context-link population");
		}
	}

	private String sourceStatus(final Source source){
		if(isUnassigned(source)){return "UNASSIGNED_SOURCE";}
		if(!source.directScored && size(source.masks)>0){return "MASKED_NO_SCORED_PARENT";}
		if(!source.directScored && size(source.diagnosticScored)>0){
			return "ONLY_WEAK_OR_BRIDGE_SCORED_PATH";
		}
		if(source.shape==PAIR){
			final Parent parent=parents.get(source.parents.get(0));
			final int kind=source.parentKinds.get(0);
			final boolean same=source.start==parent.start && source.stop==parent.stop;
			return PAIR_RESULTS[kind][same ? 0 : 1];
		}
		if(source.shape!=0){return SHAPES[source.shape];}
		if(size(source.masks)>0){return "MASKED_NO_SCORED_PARENT";}
		final int rejected=(1<<6)|(1<<8)|(1<<9)|(1<<10);
		if((source.dispositions&~rejected)==0){return "REJECTED_NO_SCORED_PARENT";}
		if(source.dispositions==(1<<1)){return "UNASSESSABLE_SOURCE";}
		if((source.dispositions&((1<<2)|(1<<7)))!=0){return "PRIOR_SOURCE_WITHOUT_EXACT_PARENT_LINK";}
		return "UNRESOLVED_NO_SCORED_PARENT";
	}

	private static String parentStatus(final Parent parent){
		if(parent.mask){return "MASKED_PARENT";}
		if(parent.prior && parent.witnesses==0){return "HITLESS_PROTECTED_PRIOR";}
		if(size(parent.sources)>0){return "SOURCE_LINKED_SCORED_PARENT";}
		if(parent.prior){return "UNLINKED_PROTECTED_PRIOR";}
		if(size(parent.unassigned)>0){return "UNASSIGNED_SOURCE_LINK_ONLY";}
		if(size(parent.otherFamilies)>0){return "OTHER_FAMILY_SOURCE_LINK_ONLY";}
		if(size(parent.diagnosticSources)>0){return "WEAK_OR_BRIDGE_SOURCE_LINK_ONLY";}
		return "CM_PARENT_WITHOUT_SOURCE_LINK";
	}

	private void writeShifts(){
		try(Output out=new Output(prefix+".shifts.tsv")){
			out.writer.println("family\tend\tchanged_unique_qualifying_group_pairs\tsigned_min\tsigned_median\tsigned_p90\tsigned_max\tabsolute_median\tabsolute_p90\tabsolute_max\tdefinition");
			for(final Family family:families){
				for(int end=0; end<2; end++){
					final IntList list=end==0 ? family.shift5 : family.shift3;
					buffer.clear().append(name(family.name)).tab().append(end==0 ? "5prime" : "3prime").tab().append(list.size);
					if(list.size==0){buffer.append("\tNA\tNA\tNA\tNA\tNA\tNA\tNA");}
					else{
						Arrays.sort(list.array, 0, list.size);
						buffer.tab().append(list.get(0)).tab().append(median(list), 1).tab().append(p90(list)).tab().append(list.get(list.size-1));
						for(int i=0; i<list.size; i++){list.array[i]=Math.abs(list.array[i]);}
						Arrays.sort(list.array, 0, list.size);
						buffer.tab().append(median(list), 1).tab().append(p90(list)).tab().append(list.get(list.size-1));
					}
					out.writer.print(buffer.append("\tPARENT_MINUS_SOURCE_ORIENTED_POSITIVE_DOWNSTREAM_MEDIAN_MIDDLE_MEAN_P90_NEAREST_RANK\n"));
				}
			}
		}
	}

	/*----------------       Byte and primitive helpers       ----------------*/

	private int symbol(final int term){
		parser.setBounds(term);
		lookup.set(parser.line(), parser.a(), parser.b());
		return intern();
	}
	private int symbol(final String value){
		BuildNcrnaReference.field(encoded.clear(), value);
		lookup.set(encoded.array, 0, encoded.length());
		return intern();
	}
	private int intern(){
		int id=symbols.get(lookup);
		if(id<0){
			id=names.size();
			final Target retained=lookup.copy();
			symbols.put(retained, id);
			names.add(retained.bytes);
		}
		return id;
	}
	private byte[] name(final int id){return names.get(id);}
	private Target probe(final ByteBuilder bytes){
		lookup.set(bytes.array, 0, bytes.length());
		return lookup;
	}
	private int input(final int term){
		final int id=symbol(term);
		require(inputCounts.contains(id), "Saved metadata input is absent from the reference manifest");
		return id;
	}
	private Family family(final int id){
		int index=familyIds.get(id);
		if(index<0){
			index=families.size();
			familyIds.put(id, index);
			families.add(new Family(id));
		}
		return families.get(index);
	}
	private Event event(final int id){
		final int index=eventIds.get(id);
		if(index<0){throw new IllegalArgumentException("Saved link names an absent event: "+id);}
		return events.get(index);
	}
	private Input open(final String name, final String header){
		final String path=new File(directory, name).getPath();
		final Input in=new Input(path);
		try{RrnaResourceIO.header(in.reader, header, path);}
		catch(RuntimeException | Error failure){
			try{in.close();}catch(RuntimeException | Error close){failure.addSuppressed(close);}
			throw failure;
		}
		return in;
	}
	private void parse(final byte[] line, final int terms){
		parser.set(line);
		require(parser.terms()==terms, "Unexpected saved comparison schema; required fields="+terms);
	}
	private int positive(final int term){
		final int value=parser.parseInt(term);
		require(value>0, "Saved row/event IDs must be positive");
		return value;
	}
	private char strand(final int term){
		parser.setBounds(term);
		require(parser.b()-parser.a()==1, "Strand must be one character");
		final char value=(char)parser.line()[parser.a()];
		require(value=='+' || value=='-' || value=='.' || value=='?', "Unknown saved strand");
		return value;
	}
	private boolean bool(final int term){
		require(parser.termEquals("true", term) || parser.termEquals("false", term), "Saved boolean must be explicit true/false");
		return parser.termEquals("true", term);
	}
	private int state(final String[] states, final int term){
		for(int i=0; i<states.length; i++){if(parser.termEquals(states[i], term)){return i;}}
		throw new IllegalArgumentException("Unknown saved state in field "+term);
	}
	private static int state(final String[] states, final String value){
		for(int i=0; i<states.length; i++){if(states[i].equals(value)){return i;}}
		throw new IllegalArgumentException("Unknown output state: "+value);
	}
	private ByteBuilder appendTerm(final ByteBuilder out, final int term){
		parser.setBounds(term);
		return out.append(parser.line(), parser.a(), parser.b()-parser.a());
	}
	private ByteBuilder geometry(final ByteBuilder out, final Node node){
		return out.append(name(node.input)).tab().append(name(node.target)).tab().append(name(node.family)).tab()
			.append(node.start).tab().append(node.stop).tab().append(node.strand);
	}
	private boolean isUnassigned(final Source source){return Arrays.equals(name(source.family), UNASSIGNED);}
	private static IntList add(IntList list, final int value){
		if(list==null){list=new IntList(1);}
		if(!contains(list, value)){list.add(value);}
		return list;
	}
	/** Edge kinds survive duplicate source/model paths. Mixed direct kinds remain
	 * explicitly separate rather than treating prior bounds as CM-group bounds. */
	private static void addDirect(final Source source, final int parent, final int kind){
		assert(kind==PRIOR_KIND || kind==GROUP_KIND) : "Only allow-listed direct associations carry a scored graph edge";
		if(source.parents==null){
			source.parents=new IntList(1);
			source.parentKinds=new IntList(1);
		}
		for(int i=0; i<source.parents.size; i++){
			if(source.parents.get(i)==parent){
				source.parentKinds.array[i]|=kind;
				return;
			}
		}
		source.parents.add(parent);
		source.parentKinds.add(kind);
	}
	private static boolean contains(final IntList list, final int value){
		if(list!=null){for(int i=0; i<list.size; i++){if(list.get(i)==value){return true;}}}
		return false;
	}
	private static int size(final IntList list){return list==null ? 0 : list.size;}
	private static long sum(final long[] values){
		long sum=0;
		for(final long value:values){sum=Math.addExact(sum, value);}
		return sum;
	}
	private static double median(final IntList list){
		assert(list.size>0) : "Median has no value for an empty changed-pair population";
		return ((long)list.get((list.size-1)/2)+list.get(list.size/2))*.5;
	}
	private static int p90(final IntList list){return list.get((int)((9L*list.size+9)/10)-1);}
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private static class Node{
		Node(final int input_, final int target_, final int family_, final int start_, final int stop_, final char strand_){
			input=input_;
			target=target_;
			family=family_;
			start=start_;
			stop=stop_;
			strand=strand_;
		}
		final int input, target, family, start, stop;
		final char strand;
		int component, shape;
	}
	private static final class Source extends Node{
		Source(final int input, final int target, final int family, final int start, final int stop, final char strand){
			super(input, target, family, start, stop, strand);
		}
		int aliases, dispositions;
		boolean directScored;
		long eventLinks, unparentedEvents;
		IntList parents, parentKinds, masks, unassignedParents, otherParents, diagnosticScored;
	}
	private static final class Parent extends Node{
		Parent(final int input, final int target, final int family, final int start, final int stop, final char strand, final int name_, final boolean mask_){
			super(input, target, family, start, stop, strand);
			name=name_;
			mask=mask_;
		}
		final int name;
		final boolean mask;
		boolean prior;
		int witnesses, priorAssociations;
		IntList sources, unassigned, otherFamilies, diagnosticSources;
	}
	private static final class Event extends Node{
		Event(final int input, final int target, final int family, final int start, final int stop, final char strand,
			final int model_, final int table_, final long line_, final int tier_){
			super(input, target, family, start, stop, strand);
			model=model_;
			table=table_;
			line=line_;
			tier=tier_;
		}
		final int model, table, tier;
		final long line;
		IntList parents, associations;
	}
	private static final class Family{
		Family(final int name_){name=name_;}
		final int name;
		boolean checked, sourceChecked;
		long sources;
		long contextLinks;
		final long[] observed=new long[11];
		final long[] sourceCounts=new long[SOURCE_STATES.length];
		final long[] statuses=new long[SOURCE_RESULTS.length];
		final long[] parentStatuses=new long[PARENT_RESULTS.length];
		final long[] lengths=new long[3];
		final long[] diagnosticPaths=new long[2];
		final IntList shift5=new IntList();
		final IntList shift3=new IntList();
	}
	private static final class Population{
		Population(final int input_){input=input_;}
		final int input;
		final long[] counts=new long[3];
		boolean checked;
	}

	private final String directory, prefix;
	private final LineParser1 parser=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder();
	private final ByteBuilder key=new ByteBuilder();
	private final ByteBuilder encoded=new ByteBuilder();
	private final Target lookup=new Target();
	private final ObjectIntMap<Target> symbols=new ObjectIntMap<Target>(Target.class);
	private final ObjectIntMap<Target> parentKeys=new ObjectIntMap<Target>(Target.class);
	private final ObjectIntMap<Target> sourceKeys=new ObjectIntMap<Target>(Target.class);
	private final ObjectIntMap<Target> fileKeys=new ObjectIntMap<Target>(Target.class);
	private final IntHashMap eventIds=new IntHashMap();
	private final IntHashMap inputCounts=new IntHashMap();
	private final IntHashMap familyIds=new IntHashMap();
	private final ArrayList<byte[]> names=new ArrayList<byte[]>();
	private final ArrayList<Parent> parents=new ArrayList<Parent>();
	private final ArrayList<Source> sources=new ArrayList<Source>();
	private final ArrayList<Event> events=new ArrayList<Event>();
	private final ArrayList<Family> families=new ArrayList<Family>();
	private final ArrayList<Population> inputs=new ArrayList<Population>();
	private final ArrayList<LongIntMap> files=new ArrayList<LongIntMap>();
	private int components;
	private long contextLinks;
	private static final int PAIR=1, SPLIT=2, MERGE=3, COMPLEX=4;
	private static final int PRIOR_KIND=1, GROUP_KIND=2;
	private static final String[] PAIR_KINDS={"NONE", "PROTECTED_PRIOR_BOUNDS", "QUALIFYING_GROUP", "MIXED_DIRECT_KINDS"};
	private static final String[][] PAIR_RESULTS={null,
		{"SAME_BOUNDS_PROTECTED_PRIOR", "CHANGED_BOUNDS_PROTECTED_PRIOR"},
		{"SAME_BOUNDS_QUALIFYING_GROUP", "CHANGED_BOUNDS_QUALIFYING_GROUP"},
		{"SAME_BOUNDS_MIXED_DIRECT", "CHANGED_BOUNDS_MIXED_DIRECT"}};
	private static final byte[] UNASSIGNED="UNASSIGNED".getBytes(StandardCharsets.US_ASCII);
	private static final String[] SHAPES={"NONE", "UNIQUE_PAIR", "SPLIT", "MERGE", "COMPLEX"};
	private static final String[] AGREEMENTS={"SAME_FAMILY", "OTHER_FAMILY", "UNASSIGNED_SOURCE"};
	private static final String[] ASSOCIATIONS={"SAME_FAMILY_OVERLAP_NOT_RESCORING_PRIOR_BOUNDS", "TANDEM_BRIDGE_WITNESS_NOT_REPRESENTATIVE", "WEAK_OVERLAP_WITNESS_NOT_GROUPING_SUPPORT", "QUALIFYING_SAME_FAMILY_GROUP_MEMBER"};
	private static final String[] SOURCE_RESULTS={"SAME_BOUNDS_QUALIFYING_GROUP", "CHANGED_BOUNDS_QUALIFYING_GROUP", "SAME_BOUNDS_PROTECTED_PRIOR", "CHANGED_BOUNDS_PROTECTED_PRIOR", "SAME_BOUNDS_MIXED_DIRECT", "CHANGED_BOUNDS_MIXED_DIRECT", "SPLIT", "MERGE", "COMPLEX", "MASKED_NO_SCORED_PARENT", "REJECTED_NO_SCORED_PARENT", "UNASSIGNED_SOURCE", "UNASSESSABLE_SOURCE", "PRIOR_SOURCE_WITHOUT_EXACT_PARENT_LINK", "UNRESOLVED_NO_SCORED_PARENT", "ONLY_WEAK_OR_BRIDGE_SCORED_PATH"};
	private static final String[] PARENT_RESULTS={"SOURCE_LINKED_SCORED_PARENT", "CM_PARENT_WITHOUT_SOURCE_LINK", "UNASSIGNED_SOURCE_LINK_ONLY", "OTHER_FAMILY_SOURCE_LINK_ONLY", "HITLESS_PROTECTED_PRIOR", "UNLINKED_PROTECTED_PRIOR", "MASKED_PARENT", "WEAK_OR_BRIDGE_SOURCE_LINK_ONLY"};
	private static final String[] EVENT_STATES={"CONFIDENT", "MARGINAL", "BELOW_MARGINAL", "UNCERTAIN_FLOOR", "REJECTED_WRONG_FAMILY", "COVERED_BY_PROTECTED_PRIOR", "INVALID_INTERVAL"};
	private static final String[] SOURCE_STATES={"CONTEXT_PROVENANCE", "UNASSESSABLE_SOURCE", "PROTECTED_PRIOR_SOURCE", "SUPPORTED_CONFIDENT_EVENT", "SUPPORTED_MARGINAL_EVENT", "UNCERTAINTY_UNMASKED", "REMOVED_REVIEWED_FAMILY_CONFLICT", "PROVENANCE_NEAR_PROTECTED_PRIOR", "REMOVED_BELOW_MARGINAL", "REMOVED_NO_CM_SUPPORT", "REMOVED_OTHER_FAMILY_ONLY"};
	private static final String SOURCE_HEADER="input_id\tgff\tsource_line\ttarget\tsource_type\tstart\tstop\tstrand\tdeclared_family\tselection_bits\tdisposition\tsource_flags\tqualification";
	private static final String LINK_HEADER="input_id\tgff\tsource_line\tevent_id\tsame_strand_valid\tfamily_agreement\tevent_tier";
	private static final String WITNESS_HEADER="annotation_id\tinput_id\tevidence_kind\tevent_id\tnative_table\tnative_line\tmodel\tstart\tstop\tprior_line\tassociation";
	private static final String EVENT_HEADER="event_id\tinput_id\ttarget\tstart\tstop\tstrand\tmodel\tfamily\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tleft_edge\tright_edge\tprinted_bits\tGA_bits\tpercent_GA\tnative_included\tprinted_GA_disagreement\tmarginal_percent\tprinted_floor_pass\tfloor_uncertainty\ttier\treason\tsource_flags\tdecision_flags\tnative_table\tnative_line\tscore_basis\tpercent_precision";
	private static final String FAMILY_HEADER="family\traw_events\tCONFIDENT\tMARGINAL\tBELOW_MARGINAL\tUNCERTAIN_FLOOR\tREJECTED_WRONG_FAMILY\tCOVERED_BY_PROTECTED_PRIOR\tINVALID_INTERVAL\tscored_annotations\tmasked_annotations\tprotected_priors";
	private static final String SOURCE_COUNT_HEADER="declared_family\tsource_rows\t"+String.join("\t", SOURCE_STATES);
	private static final String INPUT_HEADER="input_id\tscored_annotations\tmasked_annotations\tprotected_priors\tstatus";
	static final String[] INPUTS={"source_provenance.tsv.gz", "source_event_links.tsv.gz", "event_provenance.tsv.gz", "annotation_witnesses.tsv.gz", "family_counts.tsv", "source_family_counts.tsv", "summary.tsv"};
	static final String[] SUFFIXES={".source_rows.tsv.gz", ".paths.tsv.gz", ".components.tsv.gz", ".sources.tsv.gz", ".parents.tsv.gz", ".shifts.tsv", ".summary.tsv"};
}
