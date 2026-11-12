package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.PriorityQueue;

import map.ObjectIntMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import prok.GradeNcrnaReference.Binding;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaEvidenceCatalog.Membership;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Selects original input IDs from a completed reference inventory. No sequence,
 * GFF payload, caller output or native scan is opened. Reference-member support
 * is a recall stratum, not a query-subtype precision label. Input identities are
 * never collapsed by taxonomy or a guessed equivalence. Hash authentication is
 * the project recipe's responsibility, independently of these consistency checks.
 * Extended policies bound total original IDs and fill declared member quotas
 * before the ordinary balanced positive fill. Old policies keep their behavior.
 * @author Ganyu, Raiden
 */
final class NcrnaPanelSelector{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	static void process(final String census, final String memberships, final String reference,
		final String inventory, final String members, final String policy, final String taxonomy,
		final String family, final String prefix){
		final NcrnaPanelSelector selector=new NcrnaPanelSelector(census, memberships, reference,
			inventory, members, policy, taxonomy, family, prefix);
		selector.run();
	}

	private NcrnaPanelSelector(final String census_, final String memberships_, final String reference_,
		final String inventory_, final String members_, final String policy_, final String taxonomy_,
		final String family_, final String prefix_){
		census=census_;
		memberships=memberships_;
		reference=reference_;
		inventory=inventory_;
		members=members_;
		policy=policy_;
		taxonomy=taxonomy_;
		family=family_;
		prefix=prefix_;
		final NcrnaEvidenceCatalog catalog=NcrnaEvidenceCatalog.load(census, memberships);
		for(int i=0; i<catalog.size(); i++){
			final Candidate candidate=new Candidate(catalog.get(i));
			byId.put(candidate.genome.id, candidate);
			candidates.add(candidate);
		}
		Collections.sort(candidates, ID_ORDER);
	}

	private void run(){
		readMembers();
		readPolicy();
		readProvenance();
		readInputs();
		readSupport();
		readReference();
		readTaxonomy();
		select();
		writeResults();
		writeCompletion();
	}

	private void readMembers(){
		try(Input in=new Input(members)){
			RrnaResourceIO.header(in.reader, NcrnaSubtypeEvidence.MEMBER_HEADER, members);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				row(line, 3);
				final String f=text(0), name=text(1), model=text(2);
				final String key=f+"\t"+name;
				require(!memberIndex.contains(key) && !models.contains(model), "Family/member pairs and reference models must be unique");
				require(memberList.size()<30, "Member conservation uses at most thirty positive bits");
				if(!familyIndex.contains(f)){
					familyIndex.put(f, families.size());
					families.add(f);
					familyTokens.add(new Token(f));
				}
				memberIndex.put(key, memberList.size());
				models.add(model);
				memberList.add(new Member(f, name, model));
			}
		}
		require(!memberList.isEmpty() && familyIndex.contains(family), "Selected family needs an explicit complete inventory member map");
		for(final Candidate c:candidates){c.memberScored=new long[memberList.size()];}
	}

	private void readPolicy(){
		boolean found=false;
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try(Input in=new Input(policy)){
			final byte[] first=in.reader.nextLine();
			require(first!=null, "Missing panel policy header");
			final String header=new String(first, StandardCharsets.UTF_8);
			final boolean extended=TOTAL_POLICY_HEADER.equals(header);
			require(extended || POLICY_HEADER.equals(header), "Unknown panel policy header");
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				row(line, extended ? 9 : 7);
				final String f=text(0);
				require(familyIndex.contains(f) && seen.add(f), "Policy families must be unique declared inventory families");
				final int count=parser.parseInt(1), cap=parser.parseInt(2);
				final int total=extended ? Integer.parseInt(text(7)) : 0;
				require(total>=0 && (total>0 ? count==0 && cap>=total : count>0 && cap>=count),
					"Use a positive target within its ceiling, or target zero with an explicit bounded total");
				final boolean prior=bool(5), multi=bool(6);
				final String[] protectedNames=items(text(3)), sources=items(text(4));
				int bits=0;
				for(final String name:protectedNames){
					final String key=f+"\t"+name;
					require(memberIndex.contains(key), "Protected member is absent from its inventory family: "+key);
					final int index=memberIndex.get(key);
					require(memberList.get(index).family.equals(f) && (bits&(1<<index))==0,
						"Protected member must belong once to its declared family: "+name);
					bits|=1<<index;
				}
				if(!f.equals(family)){continue;}
				found=true;
				target=count;
				ceiling=cap;
				protectPrior=prior;
				protectMulti=multi;
				protectedBits=bits;
				totalInputs=total;
				for(final String source:sources){require(retainedSources.add(source), "Repeated retained source: "+source);}
				if(totalInputs>0){
					require(bits==0 && !prior && !multi, "Exact-total coarse policies retain fixed controls; full minority/prior protection belongs in validation");
					readMinimums(text(8));
				}else{require(!extended || parser.termEquals(".", 8), "Member minima require an explicit total-input budget");}
			}
		}
		require(found && seen.size()==families.size(), "Policy must cover every mapped family exactly once");
	}

	/** The comma-list order is the declared, deterministic quota-fill order. */
	private void readMinimums(final String specification){
		final String[] entries=items(specification);
		minimumInputs=new int[memberList.size()];
		quotaPosition=new int[memberList.size()];
		quotaOrder=new int[entries.length];
		int expected=0;
		for(final Member member:memberList){if(member.family.equals(family)){expected++;}}
		require(entries.length==expected, "Exact-total selection requires a minimum for every family member");
		for(int i=0; i<entries.length; i++){
			final String entry=entries[i];
			final int separator=entry.indexOf(':');
			require(separator>0 && separator==entry.lastIndexOf(':'), "Member minima use member:input_count");
			final String name=entry.substring(0, separator);
			final String key=family+"\t"+name;
			require(memberIndex.contains(key), "Unknown quota member: family="+family+", member="+name);
			final int index=memberIndex.get(key);
			require(quotaPosition[index]==0, "Repeated quota member: family="+family+", member="+name);
			final int count=Integer.parseInt(entry.substring(separator+1));
			require(count>0 && count<=totalInputs, "Member minimum must fit the explicit total-input budget: "+name);
			minimumInputs[index]=count;
			quotaPosition[index]=i+1;
			quotaOrder[i]=index;
		}
	}

	/** The inventory records its original metadata bindings and completion last. */
	private void readProvenance(){
		final String path=inventory+".provenance.tsv";
		final ObjectMap<String, String> values=new ObjectMap<String, String>(String.class, String.class);
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, "role\tpath_or_value", path);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				row(line, 2);
				final String key=text(0);
				require(!values.contains(key), "Repeated inventory provenance role: "+key);
				values.put(key, decoded(1));
			}
		}
		require(census.equals(values.get("census")) && memberships.equals(values.get("memberships"))
			&& reference.equals(values.get("reference")) && members.equals(values.get("members")),
			"Inventory must bind these exact original catalog/reference/member inputs");
		require("REFERENCE_INVENTORY_COMPLETE_NO_SELECTION_OR_ADJUDICATION".equals(values.get("status")),
			"Selection requires a completed inventory, never a partial failed run");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Inventory Counts       ----------------*/
	/*--------------------------------------------------------------*/

	private void readInputs(){
		final String path=inventory+".inputs.tsv";
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaReferenceInventory.INPUT_HEADER, path);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				row(line, 9);
				final Candidate c=candidate(0);
				final int f=inventoryFamily(1), bit=1<<f;
				require((c.familySeen&bit)==0, "Repeated input/family inventory row", c.genome.id);
				c.familySeen|=bit;
				final long scored=nonnegative(2), masked=nonnegative(3), prior=nonnegative(7), multi=nonnegative(8);
				require(bool(4)==(scored>0) && bool(5)==(scored==0 && masked>0) && bool(6)==(scored==0 && masked==0)
					&& prior<=scored && multi<=scored, "Inventory flags/counts must agree with original parent counts");
				if(!families.get(f).equals(family)){continue;}
				c.scored=scored;
				c.masked=masked;
				if(protectPrior && prior>0){c.reasons|=PRIOR;}
				if(protectMulti && multi>0){c.reasons|=MULTI;}
			}
		}
		final int expected=(1<<families.size())-1;
		for(final Candidate c:candidates){require(c.familySeen==expected, "Every original input needs every inventory family", c.genome.id);}
	}

	private void readSupport(){
		final String path=inventory+".members.tsv";
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaReferenceInventory.MEMBER_HEADER, path);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				row(line, 11);
				final Candidate c=candidate(0);
				final int f=inventoryFamily(1), m=inventoryMember(f, 2), bit=1<<m;
				final Member member=memberList.get(m);
				require(inventoryEquals(member.modelToken, 3) && (c.memberSeen&bit)==0,
					"Repeated or misbound input/member row", c.genome.id);
				c.memberSeen|=bit;
				final long scored=nonnegative(4);
				nonnegative(5);
				require(bool(6)==(scored>0), "Qualified member-positive flag must agree with its scored support");
				c.memberScored[m]=scored;
				if(member.family.equals(family)){
					// NcrnaSubtypeEvidence.support forbids confident support on masks;
					// each member therefore covers a subset of the scored union parents.
					require(scored<=c.scored, "Qualified member support exceeds the scored family union", c.genome.id);
					c.scoredSupport=Math.addExact(c.scoredSupport, scored);
					if(scored>0 && (protectedBits&bit)!=0){c.reasons|=MINORITY;}
				}
			}
		}
		final int expected=(1<<memberList.size())-1;
		for(final Candidate c:candidates){
			require(c.memberSeen==expected, "Every original input needs every declared inventory member", c.genome.id);
			// Every scored union parent has at least one qualified member in the
			// complete map. This is a lower bound, not equality: support overlaps.
			require(c.scoredSupport>=c.scored, "Scored family union lacks sufficient qualified member support", c.genome.id);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Original Bindings       ----------------*/
	/*--------------------------------------------------------------*/

	private void readReference(){
		int count=0;
		for(final Binding binding:GradeNcrnaReference.loadManifest(reference)){
			final Candidate c=byId.get(binding.id);
			require(c!=null && c.gff==null && c.genome.fasta.equals(binding.fasta), "Reference must bind each original input exactly once");
			c.gff=binding.gff;
			count++;
		}
		require(count==candidates.size(), "Selection consumes the complete registered reference population");
	}

	private void readTaxonomy(){
		final ObjectMap<String, Taxon> paths=new ObjectMap<String, Taxon>(String.class, Taxon.class);
		final ObjectSet<String> maps=new ObjectSet<String>(String.class);
		try(Input in=new Input(taxonomy)){
			RrnaResourceIO.header(in.reader, TAXONOMY_HEADER, taxonomy);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				row(line, 4);
				final String schema=text(0), path=text(1), root=text(2), hash=text(3);
				require(new File(path).isAbsolute() && new File(root).isAbsolute() && maps.add(path)
					&& SeedOffsetTable.validSha80(hash), "Taxonomy inputs need unique absolute paths/roots and explicit sha80 identities");
				readTaxonomyMap(schema, path, root, paths);
			}
		}
		require(!maps.isEmpty(), "Supply explicit taxonomy sources; unknown joins remain visible");
		final ObjectSet<String> foundSources=new ObjectSet<String>(String.class);
		for(final Candidate c:candidates){
			Taxon taxon=paths.get(c.genome.fasta);
			String taxid=null;
			boolean conflict=false;
			for(int i=0; i<c.genome.memberships(); i++){
				final Membership m=c.genome.membership(i);
				if(retainedSources.contains(m.source)){
					c.reasons|=CONTROL;
					foundSources.add(m.source);
				}
				if(!m.taxid.equals("NA")){
					if(taxid==null){taxid=m.taxid;}
					else if(!taxid.equals(m.taxid)){conflict=true;}
				}
				final Taxon alias=paths.get(m.suppliedFasta);
				if(alias!=null){
					if(taxon==null){taxon=alias;}
					else if(!taxon.same(alias)){conflict=true;}
				}
			}
			c.taxonomyConflict=conflict;
			c.taxon=taxon==null || conflict ? UNKNOWN : taxon;
			c.taxid=taxid==null || conflict ? "NA" : taxid;
			// Unknown IDs remain distinct inputs; this does not assert species independence.
			c.balanceKey=taxid==null ? "input:"+c.genome.id : "taxid:"+taxid;
			require(!conflict || c.scored==0 && c.reasons==0,
				"Conflicting taxonomy/source TaxID claims need resolution before selection", c.genome.id);
		}
		require(foundSources.size()==retainedSources.size(), "Every requested retained source must exist in the complete catalog");
	}

	private void readTaxonomyMap(final String schema, final String path, final String root,
		final ObjectMap<String, Taxon> paths){
		final boolean gtdb=schema.equals("gtdb_domain_phylum_basename_v1");
		require(gtdb || schema.equals("x4_domain_phylum_basename_source_v1"), "Unknown taxonomy schema: "+schema);
		// A separate parser preserves the outer manifest reader's fields.
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==(gtdb ? 3 : 4), "Taxonomy source row has the wrong declared schema", path);
				final String domain=string(p, 0), phylum=string(p, 1), basename=string(p, 2);
				final String source=gtdb ? domain : string(p, 3);
				require(domain.equals("bacteria") || domain.equals("archaea"), "Taxonomy map domain must be explicit bacteria or archaea");
				require(!basename.contains("/") && !basename.equals(".") && !basename.equals("..")
					&& (gtdb || source.equals("bacteria4") || source.equals("archaea4")), "Taxonomy basename/source cannot escape its declared root");
				final String fasta=root+"/"+source+"/"+basename+".fna.gz";
				final Taxon value=new Taxon(domain, phylum);
				final Taxon previous=paths.get(fasta);
				require(previous==null || previous.same(value), "Conflicting exact-source taxonomy", fasta);
				if(previous==null){paths.put(fasta, value);}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Selection           ----------------*/
	/*--------------------------------------------------------------*/

	/** Protect whole IDs first, then round-robin phyla and least-represented TaxIDs.
	 * Queues reuse one node per taxid group: no object allocation per selection. */
	private void select(){
		for(final Candidate c:candidates){
			c.selected=c.reasons!=0;
			if((c.reasons&CONTROL)!=0){
				fixedControls++;
				if(c.scored>0){positiveControls++;}else{nonpositiveControls++;}
			}
			if(c.scored==0){continue;}
			eligiblePositive++;
			if(c.selected){selectedPositive++;}
		}
		if(totalInputs>0){
			require(fixedControls<=totalInputs, "Fixed controls exceed the exact total-input budget: family="+family
				+", required_control_inputs="+fixedControls+", total_input_budget="+totalInputs);
			target=totalInputs-nonpositiveControls;
			require(target>0 && selectedPositive<=target, "Mandatory positives exceed the derived exact positive budget: family="+family
				+", required_positive_inputs="+selectedPositive+", positive_budget="+target+", nonpositive_controls="+nonpositiveControls+", total_input_budget="+totalInputs);
		}
		require(selectedPositive<=ceiling, "Mandatory positive inputs exceed the explicit ceiling; do not drop protection: "+selectedPositive+"/"+ceiling);
		final int wanted=Math.max(target, selectedPositive);
		require(eligiblePositive>=wanted, "Insufficient positive original inputs for the declared target: "+eligiblePositive+"/"+wanted);
		for(final int member:quotaOrder){fill(member, minimumInputs[member], wanted);}
		fill(-1, wanted, wanted);
	}

	/** Rebuild queues after each quota so newly credited overlapping IDs never
	 * remain in another queue's available list. All unchosen inputs remain usable.
	 */
	private void fill(final int member, final int wanted, final int positiveBudget){
		final ObjectMap<String, Group> byLineage=new ObjectMap<String, Group>(String.class, Group.class);
		final ArrayList<Group> groups=new ArrayList<Group>();
		int current=0, eligible=0;
		for(final Candidate c:candidates){
			if(c.scored==0 || member>=0 && c.memberScored[member]==0){continue;}
			eligible++;
			if(c.selected){current++;}
			final String key=c.taxon.label;
			Group group=byLineage.get(key);
			if(group==null){
				group=new Group(key);
				byLineage.put(key, group);
				groups.add(group);
			}
			group.add(c);
		}
		final String name=member<0 ? family : memberList.get(member).name;
		require(eligible>=wanted, "Insufficient supported inputs for member minimum: member="+name
			+", eligible_inputs="+eligible+", required_inputs="+wanted);
		Collections.sort(groups, Comparator.comparing(g -> g.key));
		for(final Group g:groups){g.initialize();}
		while(current<wanted){
			final int before=current;
			for(final Group group:groups){
				final Candidate c=group.take();
				if(c==null){continue;}
				if(selectedPositive>=positiveBudget){
					throw new IllegalArgumentException("Member quotas cannot fit the exact positive budget: member="+name
						+", supported_inputs="+current+", required_inputs="+wanted+", selected_positive_inputs="+selectedPositive
						+", positive_budget="+positiveBudget+", remaining_slots="+(positiveBudget-selectedPositive));
				}
				assert(!c.selected && c.scored>0) : "Balanced queues contain only unselected positive original inputs";
				c.selected=true;
				c.reasons|=member<0 ? FILL : QUOTA;
				selectedPositive++;
				current++;
				if(current==wanted){break;}
			}
			require(current>before, "Eligible queues exhausted before the declared positive target");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void writeResults(){
		int selected=0, maskedOnly=0, zero=0, unknown=0;
		final long[] eligibleMembers=new long[memberList.size()], selectedMembers=new long[memberList.size()];
		final long[] eligibleLoci=new long[memberList.size()], selectedLoci=new long[memberList.size()];
		final ObjectMap<String, Composition> sourceComposition=new ObjectMap<String, Composition>(String.class, Composition.class);
		final ObjectMap<String, Composition> lineageComposition=new ObjectMap<String, Composition>(String.class, Composition.class);
		final ArrayList<Composition> compositionRows=new ArrayList<Composition>();
		// Count and validate before opening even the first panel output.
		for(final Candidate c:candidates){
			composition(lineageComposition, compositionRows, "LINEAGE", c.taxon.label, c);
			for(int i=0; i<c.genome.memberships(); i++){
				composition(sourceComposition, compositionRows, "SOURCE", c.genome.membership(i).source, c);
			}
			for(int m=0; m<memberList.size(); m++){
				if(!memberList.get(m).family.equals(family) || c.memberScored[m]==0){continue;}
				eligibleMembers[m]++;
				eligibleLoci[m]=Math.addExact(eligibleLoci[m], c.memberScored[m]);
				if(c.selected){
					selectedMembers[m]++;
					selectedLoci[m]=Math.addExact(selectedLoci[m], c.memberScored[m]);
				}
			}
			if(!c.selected){continue;}
			selected++;
			if(c.scored==0){if(c.masked>0){maskedOnly++;}else{zero++;}}
			if(c.taxon==UNKNOWN){unknown++;}
		}
		require(selected==selectedPositive+maskedOnly+zero, "Selected input partition does not conserve slots: selected="+selected
			+", positive="+selectedPositive+", masked_only="+maskedOnly+", zero_target="+zero);
		require(totalInputs==0 || selected==totalInputs, "Exact total-input budget mismatch: selected="+selected+", required="+totalInputs);
		for(int m=0; m<memberList.size(); m++){
			final Member member=memberList.get(m);
			if(!member.family.equals(family)){continue;}
			require((protectedBits&(1<<m))==0 || eligibleMembers[m]==selectedMembers[m],
				"Protected reference-member coverage mismatch: member="+member.name+", selected="+selectedMembers[m]+", required="+eligibleMembers[m]);
			require(totalInputs==0 || selectedMembers[m]>=minimumInputs[m],
				"Exact-total member minimum mismatch: member="+member.name+", selected="+selectedMembers[m]
				+", required="+(totalInputs==0 ? 0 : minimumInputs[m]));
		}
		try(Output ledger=new Output(prefix+SUFFIXES[0]); Output inputs=new Output(prefix+SUFFIXES[1]);
			Output refs=new Output(prefix+SUFFIXES[2])){
			ledger.writer.println("input_id\tfasta\tfamily\tselected\tscored_parents\tmasked_parents\treasons\tdomain\tphylum\ttaxid_claim\ttaxonomy_status\treference_gff");
			refs.writer.println(GradeNcrnaReference.MANIFEST_HEADER);
			for(final Candidate c:candidates){
				literal(buffer.clear(), c.genome.id).tab();
				literal(buffer, c.genome.fasta).tab();
				literal(buffer, family).tab().append(c.selected).tab().append(c.scored).tab().append(c.masked).tab();
				appendReasons(buffer, c).tab();
				literal(buffer, c.taxon.domain).tab();
				literal(buffer, c.taxon.phylum).tab();
				literal(buffer, c.taxid).tab();
				literal(buffer, c.taxonomyConflict ? "CONFLICT_UNSELECTED" : c.taxon==UNKNOWN ? "UNKNOWN" : "EXACT_SOURCE_JOIN").tab();
				literal(buffer, c.gff).nl();
				ledger.writer.print(buffer);
				if(!c.selected){continue;}
				literal(buffer.clear().append("in="), c.genome.fasta).nl();
				inputs.writer.print(buffer);
				literal(buffer.clear(), c.genome.id).tab();
				literal(buffer, c.genome.fasta).tab();
				literal(buffer, c.gff).nl();
				refs.writer.print(buffer);
			}
		}
		try(Output summary=new Output(prefix+SUFFIXES[3]); Output support=new Output(prefix+SUFFIXES[4])){
			summary.writer.println("family\tcatalog_inputs\teligible_positive_inputs\trequested_positive_inputs\tselected_positive_inputs\tselected_masked_only_inputs\tselected_zero_target_inputs\tselected_inputs\tselected_unknown_taxonomy\tqualification"
				+(totalInputs==0 ? "" : "\trequested_total_inputs\tfixed_control_inputs\tpositive_fixed_controls\tnonpositive_fixed_controls\tadditional_positive_inputs"));
			literal(buffer.clear(), family).tab().append(candidates.size()).tab().append(eligiblePositive).tab().append(target)
				.tab().append(selectedPositive).tab().append(maskedOnly).tab().append(zero).tab().append(selected).tab().append(unknown)
				.append("\tORIGINAL_INPUTS_NOT_PROVEN_INDEPENDENT;NO_SEQUENCE_COLLAPSE");
			if(totalInputs>0){buffer.tab().append(totalInputs).tab().append(fixedControls).tab().append(positiveControls).tab().append(nonpositiveControls).tab().append(selectedPositive-positiveControls);}
			summary.writer.print(buffer.nl());
			support.writer.println("family\treference_member\teligible_support_inputs\tselected_support_inputs\tprotected\tqualification"
				+(totalInputs==0 ? "" : "\tminimum_support_inputs\tquota_order\teligible_scored_loci\tselected_scored_loci"));
			for(int m=0; m<memberList.size(); m++){
				final Member member=memberList.get(m);
				if(!member.family.equals(family)){continue;}
				final boolean protectedMember=(protectedBits&(1<<m))!=0;
				literal(buffer.clear(), family).tab();
				literal(buffer, member.name).tab().append(eligibleMembers[m]).tab().append(selectedMembers[m]).tab().append(protectedMember)
					.append("\tNONADDITIVE_REFERENCE_SUPPORT_NOT_CALLER_SUBTYPE_PRECISION");
				if(totalInputs>0){buffer.tab().append(minimumInputs[m]).tab().append(quotaPosition[m]).tab().append(eligibleLoci[m]).tab().append(selectedLoci[m]);}
				support.writer.print(buffer.nl());
			}
		}
		Collections.sort(compositionRows, Comparator.comparing(c -> c.key));
		try(Output strata=new Output(prefix+SUFFIXES[5])){
			strata.writer.println("family\tstratum_kind\tstratum\tcatalog_inputs\teligible_positive_inputs\tselected_inputs\tselected_positive_inputs\tqualification");
			for(final Composition c:compositionRows){
				literal(buffer.clear(), family).tab();
				literal(buffer, c.kind).tab();
				literal(buffer, c.label).tab().append(c.inputs).tab().append(c.positive)
					.tab().append(c.selected).tab().append(c.selectedPositive)
					.append("\tSOURCE_MEMBERSHIPS_NONADDITIVE;ORIGINAL_INPUTS_NOT_PROVEN_INDEPENDENT\n");
				strata.writer.print(buffer);
			}
		}
	}

	private static void composition(final ObjectMap<String, Composition> map, final ArrayList<Composition> rows,
		final String kind, final String label, final Candidate candidate){
		Composition c=map.get(label);
		if(c==null){
			c=new Composition(kind, label, kind+"\t"+label);
			map.put(label, c);
			rows.add(c);
		}
		// A catalog input may repeat within a source; one input is still one vote.
		if(c.last==candidate){return;}
		c.last=candidate;
		c.inputs++;
		if(candidate.scored>0){c.positive++;}
		if(candidate.selected){
			c.selected++;
			if(candidate.scored>0){c.selectedPositive++;}
		}
		assert(c.selected<=c.inputs && c.selectedPositive<=c.positive) : "Composition counts original IDs once within each stratum";
	}

	private void writeCompletion(){
		try(Output out=new Output(prefix+SUFFIXES[6])){
			out.writer.println("role\tpath_or_value");
			for(final String[] pair:new String[][]{{"census", census}, {"memberships", memberships}, {"reference", reference},
				{"inventory", inventory}, {"members", members}, {"policy", policy}, {"taxonomy", taxonomy}, {"family", family},
				{"status", "PANEL_SELECTION_COMPLETE_METADATA_ONLY_NO_SEQUENCE_COLLAPSE"}}){
				BuildNcrnaReference.field(buffer.clear(), pair[0]).tab();
				BuildNcrnaReference.field(buffer, pair[1]).nl();
				out.writer.print(buffer);
			}
		}
	}

	private static ByteBuilder appendReasons(final ByteBuilder b, final Candidate c){
		assert(c.reasons>=0 && c.reasons<64) : "Selection reasons use only the six declared positive bits";
		if(c.reasons==0){return b.append(c.scored>0 ? "NOT_SELECTED_BALANCED_LIMIT" : "NO_SCORED_TARGET_NOT_RETAINED_CONTROL");}
		boolean first=true;
		for(int i=0; i<REASONS.length; i++){
			if((c.reasons&(1<<i))==0){continue;}
			if(!first){b.append(';');}
			b.append(REASONS[i]);
			first=false;
		}
		return b;
	}
	/*--------------------------------------------------------------*/
	/*----------------      Metadata Framing        ----------------*/
	/*--------------------------------------------------------------*/

	private void row(final byte[] line, final int fields){
		parser.set(line);
		if(parser.terms()!=fields){
			throw new IllegalArgumentException("Selection metadata has an unexpected field count: "+parser.terms()+"/"+fields);
		}
	}
	/** Native inventory rows group the original ID contiguously. Decode a changed
	 * token once, then reuse its bound object without allocating per family/member. */
	private Candidate candidate(final int column){
		if(lastCandidate!=null && tokenEquals(inputToken, column)){return lastCandidate;}
		final String id=decoded(column);
		final Candidate c=byId.get(id);
		require(c!=null, "Unknown original input ID in inventory", id);
		remember(inputToken, column);
		lastCandidate=c;
		return c;
	}

	/** Retain one decoded spelling and raw token per declared field. The native
	 * writer emits the same spelling on repeated rows. Noncanonical but valid
	 * escapes are decoded once on a cache miss, never recursively (%2525 -> %25).
	 * The input map/policy are literal TSV; only inventory/provenance is escaped. */
	private int inventoryFamily(final int column){
		for(int f=0; f<familyTokens.size(); f++){
			final Token token=familyTokens.get(f);
			if(token.seen && tokenEquals(token.raw, column)){return f;}
		}
		final String value=decoded(column);
		require(familyIndex.contains(value), "Unmapped family in inventory", value);
		final int f=familyIndex.get(value);
		remember(familyTokens.get(f), column);
		return f;
	}
	private int inventoryMember(final int f, final int column){
		final String memberFamily=families.get(f);
		for(int m=0; m<memberList.size(); m++){
			final Member member=memberList.get(m);
			if(member.family.equals(memberFamily) && member.nameToken.seen && tokenEquals(member.nameToken.raw, column)){return m;}
		}
		final String value=decoded(column);
		for(int m=0; m<memberList.size(); m++){
			final Member member=memberList.get(m);
			if(member.family.equals(memberFamily) && member.name.equals(value)){
				remember(member.nameToken, column);
				return m;
			}
		}
		throw new IllegalArgumentException("Unknown inventory family/member: "+memberFamily+"/"+value);
	}
	private boolean inventoryEquals(final Token token, final int column){
		if(token.seen && tokenEquals(token.raw, column)){return true;}
		if(!token.value.equals(decoded(column))){return false;}
		remember(token, column);
		return true;
	}
	private boolean tokenEquals(final ByteBuilder token, final int column){
		parser.setBounds(column);
		final int start=parser.a(), length=parser.b()-start;
		if(token.length()!=length){return false;}
		for(int i=0; i<length; i++){if(token.array[i]!=parser.line()[start+i]){return false;}}
		return true;
	}
	private void remember(final Token token, final int column){
		remember(token.raw, column);
		token.seen=true;
	}
	private void remember(final ByteBuilder token, final int column){
		parser.setBounds(column);
		token.clear().append(parser.line(), parser.a(), parser.b()-parser.a());
	}
	private long nonnegative(final int column){
		final long value=parser.parseLong(column);
		require(value>=0, "Inventory counts must be nonnegative");
		return value;
	}
	private boolean bool(final int column){
		require(parser.termEquals("true", column) || parser.termEquals("false", column), "Boolean inventory/policy fields must be literal true or false");
		return parser.termEquals("true", column);
	}
	private String text(final int column){return string(parser, column);}
	private static String string(final LineParser1 parser, final int column){
		parser.setBounds(column);
		final String s=new String(parser.line(), parser.a(), parser.b()-parser.a(), StandardCharsets.UTF_8);
		require(!s.isEmpty() && s.indexOf('\uFFFD')<0, "Selection fields require nonempty lossless UTF-8");
		return s;
	}
	private String decoded(final int column){
		parser.setBounds(column);
		return decoder.decode(parser.line(), parser.a(), parser.b());
	}
	private static String[] items(final String text){
		if(text.equals(".")){return new String[0];}
		final String[] values=text.split(",", -1);
		for(final String value:values){require(!value.isEmpty() && !value.equals("."), "Policy lists cannot contain empty/missing items");}
		return values;
	}
	private static ByteBuilder literal(final ByteBuilder b, final String value){
		require(value!=null && value.indexOf('\t')<0 && value.indexOf('\n')<0 && value.indexOf('\r')<0,
			"Literal config/reference fields cannot contain record separators");
		for(int i=0; i<value.length(); i++){if(value.charAt(i)>127){return b.append(value.getBytes(StandardCharsets.UTF_8));}}
		return b.append(value);
	}
	private static void require(final boolean valid, final String reason){
		if(!valid){throw new IllegalArgumentException(reason);}
	}
	private static void require(final boolean valid, final String reason, final String context){
		if(!valid){throw new IllegalArgumentException(reason+": "+context);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Retained State        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Candidate{
		Candidate(final Genome genome_){genome=genome_;}
		final Genome genome;
		String gff, taxid, balanceKey;
		Taxon taxon;
		long scored, masked, scoredSupport;
		int familySeen, memberSeen, reasons;
		boolean selected, taxonomyConflict;
		long[] memberScored;
	}
	private static final class Composition{
		Composition(final String kind_, final String label_, final String key_){
			kind=kind_;
			label=label_;
			key=key_;
		}
		final String kind, label, key;
		Candidate last;
		int inputs, positive, selected, selectedPositive;
	}
	private static final Comparator<Candidate> ID_ORDER=(a,b)->{
		final int length=Integer.compare(a.genome.id.length(), b.genome.id.length());
		return length!=0 ? length : a.genome.id.compareTo(b.genome.id);
	};
	private static final class Member{
		Member(final String f, final String n, final String m){
			family=f;
			name=n;
			model=m;
			nameToken=new Token(n);
			modelToken=new Token(m);
		}
		final String family, name, model;
		final Token nameToken, modelToken;
	}
	private static final class Token{
		Token(final String value_){value=value_;}
		final String value;
		final ByteBuilder raw=new ByteBuilder();
		boolean seen;
	}
	private static final class Taxon{
		Taxon(final String d, final String p){
			domain=d;
			phylum=p;
			label=d+"/"+p;
		}
		boolean same(final Taxon other){return domain.equals(other.domain) && phylum.equals(other.phylum);}
		final String domain, phylum, label;
	}
	private static final class Bucket{
		Candidate next(){return available.get(index);}
		final ArrayList<Candidate> available=new ArrayList<Candidate>();
		int selected, index;
	}
	private static final class Group{
		Group(final String key_){key=key_;}
		void add(final Candidate c){
			assert(c.scored>0) : "Only reference-positive original inputs enter balanced-fill queues";
			Bucket bucket=byTaxid.get(c.balanceKey);
			if(bucket==null){
				bucket=new Bucket();
				byTaxid.put(c.balanceKey, bucket);
				buckets.add(bucket);
			}
			if(c.selected){bucket.selected++;}else{bucket.available.add(c);}
		}
		void initialize(){
			assert(queue.isEmpty()) : "A lineage queue is initialized once after mandatory IDs are counted";
			for(final Bucket b:buckets){if(!b.available.isEmpty()){queue.add(b);}}
		}
		Candidate take(){
			final Bucket b=queue.poll();
			if(b==null){return null;}
			final Candidate c=b.next();
			b.index++;
			b.selected++;
			if(b.index<b.available.size()){queue.add(b);}
			return c;
		}
		final String key;
		final ObjectMap<String, Bucket> byTaxid=new ObjectMap<String, Bucket>(String.class, Bucket.class);
		final ArrayList<Bucket> buckets=new ArrayList<Bucket>();
		final PriorityQueue<Bucket> queue=new PriorityQueue<Bucket>((a,b)->{
			final int counts=Integer.compare(a.selected, b.selected);
			return counts!=0 ? counts : ID_ORDER.compare(a.next(), b.next());
		});
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final String census, memberships, reference, inventory, members, policy, taxonomy, family, prefix;
	private final ArrayList<Candidate> candidates=new ArrayList<Candidate>();
	private final ObjectMap<String, Candidate> byId=new ObjectMap<String, Candidate>(String.class, Candidate.class);
	private final ArrayList<Member> memberList=new ArrayList<Member>();
	private final ArrayList<String> families=new ArrayList<String>();
	private final ArrayList<Token> familyTokens=new ArrayList<Token>();
	private final ObjectIntMap<String> memberIndex=new ObjectIntMap<String>(String.class), familyIndex=new ObjectIntMap<String>(String.class);
	private final ObjectSet<String> models=new ObjectSet<String>(String.class), retainedSources=new ObjectSet<String>(String.class);
	private final LineParser1 parser=new LineParser1('\t');
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final ByteBuilder buffer=new ByteBuilder();
	private final ByteBuilder inputToken=new ByteBuilder();
	private Candidate lastCandidate;
	private int target, ceiling, protectedBits, eligiblePositive, selectedPositive;
	private int totalInputs, fixedControls, positiveControls, nonpositiveControls;
	private int[] quotaOrder=new int[0], minimumInputs, quotaPosition;
	private boolean protectPrior, protectMulti;
	private static final Taxon UNKNOWN=new Taxon("UNKNOWN", "UNKNOWN");
	private static final int MINORITY=1, CONTROL=2, PRIOR=4, MULTI=8, FILL=16, QUOTA=32;
	private static final String[] REASONS={"PROTECTED_REFERENCE_MEMBER", "RETAINED_SOURCE", "PROTECTED_PRIOR", "MULTI_SUPPORTED", "LINEAGE_BALANCED_FILL", "MEMBER_QUOTA_FILL"};
	static final String[] SUFFIXES={".selection.tsv", ".inputs.config", ".reference.tsv", ".summary.tsv", ".support.tsv", ".strata.tsv", ".provenance.tsv"};
	static final String POLICY_HEADER="family\tpositive_target\tpositive_ceiling\tprotected_members\tretained_sources\tkeep_protected_priors\tkeep_multi_supported";
	static final String TOTAL_POLICY_HEADER=POLICY_HEADER+"\ttotal_inputs\tmember_min_inputs";
	static final String TAXONOMY_HEADER="schema\tpath\tsource_root\texpected_sha80";
}
