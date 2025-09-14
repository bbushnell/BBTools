package sort;

import stream.Read;
import tax.TaxNode;
import tax.TaxTree;

/**
 * @author Brian Bushnell
 * @date Oct 27, 2014
 *
 */

public final class ReadComparatorTaxa extends ReadComparator {
	
	private ReadComparatorTaxa(){}
	
	@Override
	public int compare(Read r1, Read r2) {
		int x=compareInner(r1, r2);
		return ascending*x;
	}
	
	private static int compareInner(Read r1, Read r2) {
		final TaxNode a0=tree.parseNodeFromHeader(r1.id, true);
		final TaxNode b0=tree.parseNodeFromHeader(r2.id, true);
		TaxNode a=a0, b=b0;

//		if(a==null){a=tree.getNode(1);}
//		if(b==null){b=tree.getNode(1);}
		
		if(a==null || b==null){
//			System.err.println("null for "+r1.id+", "+r2.id);
			if(a==null && b==null){return ReadComparatorName.comparator.compare(r1, r2);}
			return (a==null ? 1 : -1);
		}
		
		if(a==b){
			return ReadComparatorName.comparator.compare(r1, r2);
		}
		
//		final TaxNode c=tree.commonAncestor(a, b);
//		if(c==null){
//			assert(false) : r1.id+", "+r2.id+", "+a;
//			return compareSimple(tree.highestAncestor(a), tree.highestAncestor(b));
//		}
//
//		while(a.id!=c.id && a.pid!=c.id){a=tree.getNode(a.pid);}
//		while(b.id!=c.id && c.pid!=c.id){b=tree.getNode(b.pid);}
		
		while(a.id!=a.pid && a.levelExtended<TaxTree.FAMILY_E){
			TaxNode x=tree.getNode(a.pid);
			if(x.levelExtended>TaxTree.FAMILY_E){break;}
			a=x;
		}
		while(b.id!=b.pid && b.levelExtended<TaxTree.FAMILY_E){
			TaxNode x=tree.getNode(b.pid);
			if(x.levelExtended>TaxTree.FAMILY_E){break;}
			b=x;
		}

		if(a==b){
			while(a.id!=a.pid && a.levelExtended<TaxTree.GENUS_E){
				TaxNode x=tree.getNode(a.pid);
				if(x.levelExtended>TaxTree.GENUS_E){break;}
				a=x;
			}
			while(b.id!=b.pid && b.levelExtended<TaxTree.GENUS_E){
				TaxNode x=tree.getNode(b.pid);
				if(x.levelExtended>TaxTree.GENUS_E){break;}
				b=x;
			}

			if(a==b){
				a=a0;
				b=b0;
				while(a.id!=a.pid && a.levelExtended<TaxTree.SPECIES_E){
					TaxNode x=tree.getNode(a.pid);
					if(x.levelExtended>TaxTree.SPECIES_E){break;}
					a=x;
				}
				while(b.id!=b.pid && b.levelExtended<TaxTree.SPECIES_E){
					TaxNode x=tree.getNode(b.pid);
					if(x.levelExtended>TaxTree.SPECIES_E){break;}
					b=x;
				}

				if(a==b){
					return compareSimple(a0, b0);
				}
			}
		}
		
		return compareSimple(a, b);
	}
	
	private static int compareSimple(final TaxNode a, final TaxNode b){
		if(a.levelExtended!=b.levelExtended){return a.levelExtended-b.levelExtended;}
		return a.id-b.id;
	}
		
	private int ascending=1;
	
	@Override
	public void setAscending(boolean asc){
		ascending=(asc ? 1 : -1);
	}

	public static TaxTree tree;
	
	public static final ReadComparatorTaxa comparator=new ReadComparatorTaxa();
	
}
