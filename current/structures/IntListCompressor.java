package structures;

public final class IntListCompressor {
	
	public void add(int value){
		list.add(value);
		if(list.freeSpace()==0 && lastCompression<0.75f*list.size()){
			sortAndShrink();
		}
	}
	
	public void sortAndShrink(){
		if(lastCompression>=list.size()){return;}
		list.sort();
		list.shrinkToUnique();
		lastCompression=list.size();
	}
	
	public IntList list=new IntList(4);
	private int lastCompression=0;
	
}
