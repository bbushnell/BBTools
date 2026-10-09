package tax;

import java.util.ArrayList;
import java.util.HashMap;

import fileIO.ByteFile;
import shared.Shared;
import shared.Tools;
import structures.IntHashMap;

/** Adds an auxiliary namespace of durable fake TaxIDs to a packed runtime TaxTree. */
public final class TaxTreeAux {

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Load {@code fname} and overlay its fake nodes onto {@code base}. */
	public static TaxTree apply(TaxTree base, String fname){
		if(base==null){throw new NullPointerException("base");}
		ArrayList<AuxRow> rows=load(fname);
		if(rows.isEmpty()){return base;}
		final int firstAuxIndex=base.nodes.length;
		TaxNode[] nodes=cloneNodes(base.nodes, rows.size());
		IntHashMap externalToInternal=new IntHashMap(Tools.max(2, rows.size()*2));
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
			final int internal=firstAuxIndex+i;
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

	private static TaxNode getParent(TaxTree base, TaxNode[] nodes, IntHashMap map, AuxRow row){
		TaxNode parent=getNode(base, nodes, map, row.parentTaxID);
		if(parent==null){
			throw new IllegalArgumentException("Auxiliary TaxID "+row.fakeTaxID+
					" has unknown parent "+row.parentTaxID);
		}
		return parent;
	}

	private static TaxNode getNode(TaxTree base, TaxNode[] nodes, IntHashMap map, int id){
		int internal=map.get(id);
		if(internal>=0){return nodes[internal];}
		return id<base.nodes.length ? nodes[id] : null;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class AuxRow {
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

		final int fakeTaxID, parentTaxID, level, levelExtended;
		final String name;
	}

	private static final class Header {
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

		final int fakeTaxID, parentTaxID, rank, name, minColumns;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private static int parsePositive(String s, String field){
		final long x=Long.parseLong(s);
		if(x<1 || x>Integer.MAX_VALUE){throw new IllegalArgumentException("Invalid "+field+": "+s);}
		return (int)x;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	public static final int MIN_FAKE_TAXID=2000000000;

}
