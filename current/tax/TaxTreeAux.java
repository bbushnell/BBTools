package tax;

import java.util.ArrayList;
import java.util.HashMap;

import fileIO.ByteFile;
import shared.Shared;
import shared.Tools;
import structures.IntHashMap;

/** Loads a named auxiliary TaxID table and builds a copied, packed overlay tree.
 * Fake IDs remain durable external identifiers while their nodes occupy appended
 * array slots. Existing nodes are copied; the base tree's merged-ID map is shared.
 * Parsed row order determines packed positions, not the externally visible IDs. */
public final class TaxTreeAux{

	/** Creates a stateless holder; overlay operations are static. */
	public TaxTreeAux(){}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Copies existing nodes and attaches an auxiliary hierarchy to native or auxiliary IDs.
	 * Empty tables return base unchanged. Nonempty tables validate IDs, duplicate
	 * rows, auxiliary cycles and parent presence, then propagate auxiliary rank bounds
	 * and rebuild rank indexes. Existing auxiliary IDs and their private-slot boundary
	 * are preserved in a copied map. Optional runtime name/size indexes are not copied.
	 * @param base Base topology to overlay, never null
	 * @param fname Tab-delimited auxiliary table with a named header
	 * @return Fresh overlay tree, or base for a header-only table */
	public static TaxTree apply(TaxTree base, String fname){
		if(base==null){throw new NullPointerException("base");}
		ArrayList<AuxRow> rows=load(fname);
		if(rows.isEmpty()){return base;}
		final int appendStart=base.nodes.length;
		final IntHashMap oldMap=base.externalToInternalMap;
		final int firstAuxIndex=(oldMap==null ? appendStart : base.firstAuxIndex);
		TaxNode[] nodes=cloneNodes(base.nodes, rows.size());
		final long entries=rows.size()+(oldMap==null ? 0L : oldMap.size());
		final int capacity=(int)Math.max(2L, Math.min(Integer.MAX_VALUE, 2L*entries));
		IntHashMap externalToInternal=new IntHashMap(capacity);
		if(oldMap!=null){externalToInternal.putAll(oldMap);}
		HashMap<Integer, AuxRow> rowsById=new HashMap<Integer, AuxRow>(Tools.max(2, rows.size()*2));
		for(int i=0; i<rows.size(); i++){
			AuxRow row=rows.get(i);
			if(row.fakeTaxID<MIN_FAKE_TAXID){
				throw new IllegalArgumentException("Auxiliary TaxID "+row.fakeTaxID+" is below "+MIN_FAKE_TAXID);
			}
			if(base.getNode(row.fakeTaxID, true)!=null){
				throw new IllegalArgumentException("Auxiliary TaxID "+row.fakeTaxID+" collides with the base tree");
			}
			if(rowsById.put(row.fakeTaxID, row)!=null){
				throw new IllegalArgumentException("Duplicate auxiliary TaxID "+row.fakeTaxID);
			}
			final int internal=appendStart+i;
			externalToInternal.put(row.fakeTaxID, internal);
			nodes[internal]=new TaxNode(row.fakeTaxID, row.parentTaxID, row.level, row.levelExtended, row.name);
			nodes[internal].setCanonical(TaxTree.isSimple(row.levelExtended));
		}
		validateNoAuxCycles(rows, rowsById);
		for(AuxRow row : rows){
			TaxNode parent=getParent(base, nodes, externalToInternal, row);
			if(parent.id==row.fakeTaxID){
				throw new IllegalArgumentException("Auxiliary TaxID "+row.fakeTaxID+" cannot be its own parent");
			}
			parent.numChildren++;
		}
		discussWithParents(base, nodes, externalToInternal, firstAuxIndex);
		return new TaxTree(nodes, base.mergedMap, base.minValidTaxa, base.simplify,
				base.reassign, base.skipNorank, base.inferRankLimit, Shared.threads()/2,
				externalToInternal, firstAuxIndex);
	}

	/** Reads nonblank, noncomment table rows and propagates the reader's close status.
	 * @param fname Input filename
	 * @return Parsed rows in file order; the header itself is omitted */
	private static ArrayList<AuxRow> load(String fname){
		if(fname==null || fname.length()<1){throw new IllegalArgumentException("Missing auxiliary tree file");}
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		ArrayList<AuxRow> rows=new ArrayList<AuxRow>();
		Header header=null;
		boolean error=false;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length<1 || line[0]=='#'){continue;}
				String[] split=new String(line).split("\t", -1);
				if(header==null){header=new Header(split);}
				else{rows.add(new AuxRow(split, header));}
			}
		}finally{
			error=bf.close();
		}
		if(error){throw new RuntimeException("Error reading "+fname);}
		if(header==null){throw new IllegalArgumentException("Missing header in auxiliary tree file "+fname);}
		return rows;
	}

	/** Copies existing nodes, including flags, counters and inferred bounds.
	 * Names remain shared immutable strings; extra slots remain null.
	 * @param old Source node array
	 * @param extra Number of appended auxiliary slots
	 * @return Expanded array containing independent node objects */
	private static TaxNode[] cloneNodes(TaxNode[] old, int extra){
		TaxNode[] nodes=new TaxNode[old.length+extra];
		for(int i=0; i<old.length; i++){
			TaxNode src=old[i];
			if(src!=null){
				TaxNode dst=new TaxNode(src.id, src.pid, src.level, src.levelExtended, src.name);
				dst.numChildren=src.numChildren;
				dst.minParentLevelExtended=src.minParentLevelExtended;
				dst.maxChildLevelExtended=src.maxChildLevelExtended;
				dst.countRaw=src.countRaw;
				dst.countSum=src.countSum;
				dst.setRawFlag(src.rawFlag());
				nodes[i]=dst;
			}
		}
		return nodes;
	}

	/** Rejects parent walks longer than the number of auxiliary rows.
	 * A walk stops when its parent leaves the auxiliary namespace; native topology
	 * is assumed valid and is not traversed by this check.
	 * @param rows All new auxiliary rows
	 * @param rowsById Complete mapping of new auxiliary IDs to rows */
	private static void validateNoAuxCycles(ArrayList<AuxRow> rows, HashMap<Integer, AuxRow> rowsById){
		for(AuxRow row : rows){
			AuxRow current=row;
			for(int hops=0; current!=null; hops++){
				if(hops>rows.size()){
					throw new IllegalArgumentException("Cycle detected under auxiliary TaxID "+row.fakeTaxID);
				}
				current=rowsById.get(current.parentTaxID);
			}
		}
	}

	/** Propagates rank bounds over auxiliary edges until no bound changes.
	 * @param base Base tree supplying native-ID array bounds
	 * @param nodes Copied native nodes followed by appended auxiliary nodes
	 * @param map New auxiliary IDs to packed slots
	 * @param firstAuxIndex First appended slot */
	private static void discussWithParents(TaxTree base, TaxNode[] nodes, IntHashMap map, int firstAuxIndex){
		boolean changed=true;
		while(changed){
			changed=false;
			for(int i=firstAuxIndex; i<nodes.length; i++){
				TaxNode child=nodes[i];
				TaxNode parent=getNode(base, nodes, map, child.pid);
				changed=(child.discussWithParent(parent) | changed);
			}
		}
	}

	/** Resolves a parent and rejects a missing parent with the child's durable ID.
	 * @param base Base tree supplying native-ID bounds
	 * @param nodes Working overlay nodes
	 * @param map New auxiliary-ID mapping
	 * @param row Child whose parent is required
	 * @return Resolved parent in the copied working tree */
	private static TaxNode getParent(TaxTree base, TaxNode[] nodes, IntHashMap map, AuxRow row){
		TaxNode parent=getNode(base, nodes, map, row.parentTaxID);
		if(parent==null){
			throw new IllegalArgumentException("Auxiliary TaxID "+row.fakeTaxID+
					" has unknown parent "+row.parentTaxID);
		}
		return parent;
	}

	/** Resolves a positive parent ID in the combined auxiliary map or base tree.
	 * Private packed slots are not accepted as externally supplied parent IDs.
	 * @param base Base tree resolving external IDs and hiding private slots
	 * @param nodes Working overlay nodes
	 * @param map Combined old and new auxiliary IDs to packed positions
	 * @param id Positive durable TaxID
	 * @return Resolved node, or null if the ID is absent */
	private static TaxNode getNode(TaxTree base, TaxNode[] nodes, IntHashMap map, int id){
		int internal=map.get(id);
		if(internal>=0){return nodes[internal];}
		internal=base.nodeIndex(id);
		return internal>=0 && internal<base.nodes.length ? nodes[internal] : null;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parsed auxiliary record with normal and extended rank resolved once. */
	private static final class AuxRow{
		/** Validates required fields and resolves textual or numeric rank input.
		 * Numeric ranks follow parseLevelExtended's normal-rank-number convention.
		 * @param split Complete tab fields, including empty trailing fields
		 * @param header Required-column positions */
		AuxRow(String[] split, Header header){
			if(split.length<header.minColumns){
				throw new IllegalArgumentException("Expected at least "+header.minColumns+
						" columns but found "+split.length);
			}
			fakeTaxID=parsePositive(split[header.fakeTaxID], "fake_taxid");
			parentTaxID=parsePositive(split[header.parentTaxID], "parent_taxid");
			name=split[header.name];
			if(name.length()<1){throw new IllegalArgumentException("Missing name for TaxID "+fakeTaxID);}
			String rank=split[header.rank];
			if(rank.length()<1 || (!Tools.isNumeric(rank.charAt(0)) && !TaxTree.levelMapExtendedContains(rank.toLowerCase()))){
				throw new IllegalArgumentException("Invalid rank "+rank+" for TaxID "+fakeTaxID);
			}
			levelExtended=TaxTree.parseLevelExtended(rank);
			if(levelExtended<0 || levelExtended>=TaxTree.numTaxLevelNamesExtended){
				throw new IllegalArgumentException("Invalid rank "+split[header.rank]+" for TaxID "+fakeTaxID);
			}
			level=TaxTree.extendedToLevel(levelExtended);
		}

		/** Durable child/parent IDs and resolved normal/extended ranks. */
		final int fakeTaxID, parentTaxID, level, levelExtended;
		/** Nonempty name, retained as parsed. */
		final String name;
	}

	/** Required-column positions from the first noncomment table row. */
	private static final class Header{
		/** Resolves column names case-insensitively; later duplicate names take precedence.
		 * @param split Header fields, with optional extra columns */
		Header(String[] split){
			int fakeTaxID_=-1, parentTaxID_=-1, rank_=-1, name_=-1;
			for(int i=0; i<split.length; i++){
				String s=split[i].toLowerCase();
				if(s.equals("fake_taxid")){fakeTaxID_=i;}
				else if(s.equals("parent_taxid")){parentTaxID_=i;}
				else if(s.equals("rank")){rank_=i;}
				else if(s.equals("name")){name_=i;}
			}
			if(fakeTaxID_<0 || parentTaxID_<0 || rank_<0 || name_<0){
				throw new IllegalArgumentException("Auxiliary tree header must include fake_taxid, parent_taxid, rank, and name");
			}
			fakeTaxID=fakeTaxID_;
			parentTaxID=parentTaxID_;
			rank=rank_;
			name=name_;
			minColumns=Tools.max(fakeTaxID, parentTaxID, rank, name)+1;
		}

		/** Required column positions and minimum row width covering them. */
		final int fakeTaxID, parentTaxID, rank, name, minColumns;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses a positive decimal TaxID bounded by Integer.MAX_VALUE.
	 * @param s Decimal input token
	 * @param field Field name for diagnostics
	 * @return Valid positive integer */
	private static int parsePositive(String s, String field){
		final long x=Long.parseLong(s);
		if(x<1 || x>Integer.MAX_VALUE){throw new IllegalArgumentException("Invalid "+field+": "+s);}
		return (int)x;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Lowest permitted durable auxiliary TaxID; native parent IDs may be smaller. */
	public static final int MIN_FAKE_TAXID=2000000000;

}
