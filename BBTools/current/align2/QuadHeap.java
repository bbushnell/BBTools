package align2;

public final class QuadHeap {
	
	public QuadHeap(int maxSize){
		
		int len=maxSize+1;
		if((len&1)==1){len++;} //Array size is always even.
		
		CAPACITY=maxSize;
		array=new Quad[len];
//		queue=new PriorityQueue<T>(maxSize);
	}
	
	public boolean add(Quad t){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
//		queue.add(t);
		assert(size==0 || array[size]!=null);
		size++;
		array[size]=t;
		percDown(size);
//		assert(queue.size()==size);
//		assert(queue.peek()==peek());
		//assert(testForDuplicates());
		return true;
	}
	
	public Quad peek(){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
		if(size==0){return null;}
//		assert(array[1]==queue.peek()) : size+", "+queue.size()+"\n"+
//			array[1]+"\n"+
//			array[2]+" , "+array[3]+"\n"+
//			array[4]+" , "+array[5]+" , "+array[6]+" , "+array[7]+"\n"+
//			queue.peek()+"\n";
		//assert(testForDuplicates());
		return array[1];
	}
	
	public Quad poll(){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
		if(size==0){return null;}
		Quad t=array[1];
//		assert(t==queue.poll());
		array[1]=array[size];
		array[size]=null;
		size--;
		if(size>0){percUp(1);}
//		assert(queue.size()==size);
//		assert(queue.peek()==peek());
		//assert(testForDuplicates());
		return t;
	}
	
//	private void percDownRecursive(int loc){
//		//assert(testForDuplicates());
//		assert(loc>0);
//		if(loc==1){return;}
//		int next=loc/2;
//		Quad a=array[loc];
//		Quad b=array[next];
//		assert(a!=b);
//		if(a.compareTo(b)<0){
//			array[next]=a;
//			array[loc]=b;
//			percDown(next);
//		}
//	}
//
//	private void percDown_old(int loc){
//		//assert(testForDuplicates());
//		assert(loc>0);
//
//		final Quad a=array[loc];
//
//		while(loc>1){
//			int next=loc/2;
//			Quad b=array[next];
//			assert(a!=b);
//			if(a.compareTo(b)<0){
//				array[next]=a;
//				array[loc]=b;
//				loc=next;
//			}else{return;}
//		}
//	}
	
	private void percDown(int loc){
		//assert(testForDuplicates());
		assert(loc>0);
		if(loc==1){return;}

		int next=loc/2;
		final Quad a=array[loc];
		Quad b=array[next];
		
//		while(loc>1 && (a.site<b.site || (a.site==b.site && a.column<b.column))){
		while(loc>1 && a.compareTo(b)<0){
			array[loc]=b;
			loc=next;
			next=next/2;
			b=array[next];
		}
			
		array[loc]=a;
	}
	
	private void percUp(int loc){
		//assert(testForDuplicates());
		assert(loc>0 && loc<=size) : loc+", "+size;
		int next1=loc*2;
		int next2=next1+1;
		if(next1>size){return;}
		Quad a=array[loc];
		Quad b=array[next1];
		Quad c=array[next2];
		assert(a!=b);
		assert(b!=c);
		assert(b!=null);
		//assert(testForDuplicates());
		if(c==null || b.compareTo(c)<1){
			if(a.compareTo(b)>0){
//			if((a.site>b.site || (a.site==b.site && a.column>b.column))){
				array[next1]=a;
				array[loc]=b;
				//assert(testForDuplicates());
				percUp(next1);
			}
		}else{
			if(a.compareTo(c)>0){
//			if((a.site>c.site || (a.site==c.site && a.column>c.column))){
				array[next2]=a;
				array[loc]=c;
				//assert(testForDuplicates());
				percUp(next2);
			}
		}
	}
	
	private void percUpIter(int loc){
		//assert(testForDuplicates());
		assert(loc>0 && loc<=size) : loc+", "+size;
		final Quad a=array[loc];
		//assert(testForDuplicates());

		int next1=loc*2;
		int next2=next1+1;
		
		while(next1<=size){
			
			Quad b=array[next1];
			Quad c=array[next2];
			assert(a!=b);
			assert(b!=c);
			assert(b!=null);
			
			if(c==null || b.compareTo(c)<1){
//			if(c==null || (b.site<c.site || (b.site==c.site && b.column<c.column))){
				if(a.compareTo(b)>0){
//				if((a.site>b.site || (a.site==b.site && a.column>b.column))){
//					array[next1]=a;
					array[loc]=b;
					loc=next1;
				}else{
					break;
				}
			}else{
				if(a.compareTo(c)>0){
//				if((a.site>c.site || (a.site==c.site && a.column>c.column))){
//					array[next2]=a;
					array[loc]=c;
					loc=next2;
				}else{
					break;
				}
			}
			next1=loc*2;
			next2=next1+1;
		}
		array[loc]=a;
	}
	
	public boolean isEmpty(){
//		assert((size==0) == queue.isEmpty());
		return size==0;
	}
	
	public void clear(){
//		queue.clear();
//		for(int i=1; i<=size; i++){array[i]=null;}
		size=0;
	}
	
	public int size(){
		return size;
	}
	
	public static int tier(int x){
		int leading=Integer.numberOfLeadingZeros(x);
		return 31-leading;
	}
	
	public boolean testForDuplicates(){
		for(int i=0; i<array.length; i++){
			for(int j=i+1; j<array.length; j++){
				if(array[i]!=null && array[i]==array[j]){return false;}
			}
		}
		return true;
	}
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append("[");
		for(int i=1; i<=size; i++){
			sb.append((i==1 ? "" : ", ")+array[i]);
		}
		sb.append("]");
		return sb.toString();
	}
	
	private final Quad[] array;
	private final int CAPACITY;
	private int size=0;
	
}
