#include "raylib.h"
#include <stdbool.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct Pane Pane;
typedef struct Cursor Cursor ;

struct Cursor{
	int w,h;
	double x,y,vx,vy,s;
};

struct Pane{
	Rectangle r;
	int id, p_id, pix_w, pix_h;
	Color (*c_fun)(double x,double y, Pane *p);
	Color *pix;
	Pane *parent;
	int is_drawable, is_grabable, is_movable,is_selectable;

};


Color c_fun_from_pix(double x, double y, Pane *p){
	int ix = (int) x*(p->pix_w), iy = (int) y*(p->pix_h);
	//printf("\n--------------\n(ix, iy) = (%d, %d)\n",ix,iy);
	Color c = p->pix[iy*p->pix_w + ix];
	//printf("c = (a,r,g,b) = (%d,%d,%d,%d)\n",c.a,c.r,c.g,c.b);
	return c;
}

Color c_fun_grad_rb(double x, double y,Pane *p){
	double dx = 0.5 - x, dy = 0.5 -y, d = fabs(dx) + fabs(dy); 
	return (Color){.a = 255,.r = 255*d,.g = 0,.b = 255*(1-d)};
}

Color c_fun_rand(double x, double y,Pane *p){
	double p1 = 1.0*rand()/RAND_MAX, p2 = 1.0*rand()/RAND_MAX, p3 = 1.0*rand()/RAND_MAX;
	return (Color){.a = 255,.r = 255*p1,.g = 255*p2,.b = 255*p3};
}

Color c_fun_green(double x, double y, Pane *p){
	return GREEN;
}

void draw_color(double x, double y, Color *pix, Rectangle r, int w, int h,Color c){
	int x_n = r.x + x*(r.width/w);
	int y_n = r.y + y*(r.height/h);
	pix[y_n*w + x_n] = c;
}


int min(int a, int b){
	return (a-b < 0) ? a : b;
}


int IsCursorOnPane(Pane *p,Cursor *c,int w, int h){
	int mw = c->x, mh = c->y;
	int px = p->r.x*w, py = p->r.y*h, pw = p->r.width*w, ph = p->r.height*h;

	return ((px <= mw)*(mw < px+ pw)*(py<=mh)*(mh<py+ph));
}

void IsCursorOnPanes(Pane *p,Cursor *c,int* cursor_on_panes,int np,int w, int h){
	for(int i = 0;i<np;++i){
		cursor_on_panes[i] = IsCursorOnPane(p+i,c,w,h);
	}
}

int pane_can_move(Pane *p, Cursor *c){
	double px = p->r.x, py = p->r.y, pw = p->r.width, ph = p->r.height;
	double ppx = p->parent->r.x, ppy = p->parent->r.y, ppw = p->parent->r.width, pph = p->parent->r.height;
	double n_x = px + c->vx/c->w, n_y = py +c->vy/c->h, n_w = pw + 0.0*c->s, n_h = ph + 0.0*c->s ;
	return (n_x +n_w < ppx + ppw  && n_x > ppx)*(n_y + n_h < ppy+pph && n_y > ppy );
}
void move_pane(Pane *p,Cursor *c){
	double px = p->r.x, py = p->r.y, pw = p->r.width, ph = p->r.height;
	double ppx = (p - p->id + p->p_id)->r.x, ppy = (p - p->id + p->p_id)->r.y, ppw = (p - p->id + p->p_id)->r.width, pph = (p - p->id + p->p_id)->r.height;
	double vx = c->vx/c->w*ppw, vy = c->vy/c->h*pph, vs = c->s/sqrt(c->h*c->w);
	double n_x = px + vx, n_y = py +vy, n_w = pw + vs, n_h = ph + vs ;
	if(n_x +pw <= ppx + ppw  && n_x >= ppx){
		p->r.x = n_x;
	}
	if(n_y + ph <= ppy+pph && n_y >= ppy ){
		p->r.y = n_y;}
	if(n_x + n_w < ppx+ppw && n_y + n_h < ppy+pph){
		p->r.height = n_h; p->r.width = n_w;}
}


void upd_childs(Pane *p, int np,int i,double opw, double oph, double opx, double opy){
	double px = (p+i)->r.x, py = (p+i)->r.y, pw = (p+i)->r.width, ph = (p+i)->r.height, pid = (p+i)->id;
	//double opx = p->r.x, opy = p->r.y, opw = p->r.width, oph = p->r.height;
	for(int k = 0;k<np;++k){
		Pane *chi = p+k;
		double ocx = chi->r.x, ocy = chi->r.y, ocw = chi->r.width, och = chi->r.height;
		double cx = ocx, cy = ocy, cw = ocw, ch = och;
		if(chi->p_id == pid){

			Rectangle ol_r_chi = chi->r;
			chi->r.width = (pw/opw)*cw;
			chi->r.height = (ph/oph)*ch;
			chi->r.x = (cx-opx)/opw;
			chi->r.y = (cy-opy)/oph;
			double cx = chi->r.x, cy = chi->r.y, cw = chi->r.width, ch = chi->r.height;
			chi->r.x = (1-cx)*px + cx*(px+pw);
			chi->r.y = (1-cy)*py + cy*(py+ph);
			upd_childs(p, np, k,ocw,och,ocx,ocy);
		}else{
		}
	}

	//move_pane(p+i, c);

}


double lintrp(double x, double y, double t){
	return (1-t)*x + t*y;
}

Pane *init_sub_pane(Pane *p,double x, double y, double w, double h,Color (*cf)(double x, double y, Pane *p),int n_panes,Color *pix, int pix_w, int pix_h){
	Pane *sp = malloc(sizeof(Pane));
	double spy = p->r.y*(1-y) + (p->r.y+p->r.height)*y, spx = p->r.x*(1-x) + (p->r.x+p->r.width)*x, sph = p->r.height*h, spw = p->r.width*w;
	sp->r = (Rectangle) {.width = spw,.height = sph,.y = spy,.x = spx};
	sp->parent = p;
	sp->c_fun = cf;
	sp->is_movable = 1;
	sp->is_grabable = 1;
	sp->is_selectable = 1;
	sp->is_drawable = 1;
	sp->id = n_panes;
	sp->p_id = p->id;
	sp->pix_h = pix_h;
	sp->pix_w = pix_w;
	sp->pix = pix;
	return sp;
}







void update_pane_pixels(Color *pix,int w,int h,Pane *p){
	int px = p->r.x*w, py = p->r.y*h, pw = p->r.width*w, ph = p->r.height*h;
	for(int i = py;i<min(h,py+ph);++i){
		//if(py+i > h || py +i <0){continue;};
		for(int j = px;j<min(w,px+pw);++j){
		//	if(px+j>w || px+j < 0){continue;};
			double x = 1.0*(j-px)/pw, y = 1.0*(i-py)/ph;
			pix[i*w + j] = p->c_fun(x,y,p);

		}
	}
}

void update_panes_pixels(Pane *ps, int n_panes, Color *pix,int w, int h){
	for(int i = 0;i<n_panes;++i){
		Pane *p = ps + i;
		if(p->is_drawable){
			update_pane_pixels(pix, w, h,p);
		}

	}
}



void update_pixels(Color *pixels,int w,int h){
	//some code to update the pixs;	
}



int main() {
	int w = 1000, h = 1000;
	SetConfigFlags(FLAG_WINDOW_RESIZABLE);
	InitWindow(w, h, "");
	SetTargetFPS(60);
	Image img = GenImageColor(w, h, BLACK);
	// set the image's format so it is guaranteed to be aligned with our pixel buffer format below
	//ImageFormat(&img, PIwELFORMAh_UNCOMPRESSED_R8G8B8A8);
	Texture tex = LoadTextureFromImage(img);

	Color rgba_pixels[h * w]; // 4 channels


	// make sure to set ALL ALPHA CHANNELS of the pixels (every fourth index) to 255, otherwise you'll have transparent pixels
	
	// push the rgba pixel values to the texture
	
	int a = 0;
	int scale = 1;	
	Rectangle rec = {.x = 0,.y = 0,.width = w,.height = h };
	Pane w_pane = (Pane) {.r = (Rectangle) {.x=0,.y=0,.width=1,.height=1},.is_drawable = 0,.id = 0};
	int n_panes = 2;
	Pane *p_ptr = malloc(n_panes*sizeof(Pane));
	*p_ptr = w_pane;
	w_pane.is_grabable = 1;
	w_pane.is_movable = 1;
	*(p_ptr +1) = *init_sub_pane(p_ptr, 0, 0,0.5, 0.5, c_fun_from_pix,1,&RED,1,1);
	//*(p_ptr +2) = *init_sub_pane(p_ptr+1, 0, 0,0.5, 0.5, c_fun_grad_rb,2);
	//*(p_ptr +3) = *init_sub_pane(p_ptr+1, 0.1, 0.1,0.5, 0.5, c_fun_rand,3);
	int held_idx = -1;
	int *panes_held = malloc(n_panes);
	int *panes_selected = malloc(n_panes);
	int *cursor_on_panes = malloc(n_panes);
	Color *sub_pix = &BLUE;
	Cursor cursor = {.x = 1.0*w*GetMouseX()/GetScreenWidth(),.y = 1.0*GetMouseY()/GetScreenHeight(),.vx = 0,.vy = 0,.s = 0,.h = h,.w = w};
	while(!WindowShouldClose()){
		for(int i = 0;i<h*w;++i){
			rgba_pixels[i] = WHITE;
		}
		Rectangle big_rec = {.x = 0,.y = 0,.width = GetScreenWidth(),.height = GetScreenHeight() };
		IsCursorOnPanes(p_ptr,&cursor,cursor_on_panes,n_panes,w,h);
		//printf("%d \n",cursor_on_panes[0]);
		//printf("cursor pos (x,y) = (%d,%d),held_idx = %d \n",cursor.x,cursor.y,held_idx);
		int b = IsMouseButtonDown(0);
		
		if(held_idx == -1 && b){
			for(int i = 0;i<n_panes;++i){
				if(cursor_on_panes[i] && (p_ptr+i)->is_grabable){
					held_idx = i;
				}
		};};
		if(!b){held_idx = -1;};
		//update pixels before drawin 
		// update_pixels(rgba_pixels,w,h,...)
		//p_ptr[0].r.x = 1.0*w*rand()/RAND_MAX;
		//p_ptr[0].r.y = 1.0*h*rand()/RAND_MAX;
		update_panes_pixels(p_ptr,n_panes,rgba_pixels, w, h);
		//load new pixs into texture
		UpdateTexture(tex, rgba_pixels);
		BeginDrawing();
		//reset window and draw texture
		ClearBackground(RAYWHITE);
		DrawTexturePro(tex,rec,big_rec,(Vector2){0,0},0,WHITE);
		EndDrawing();
		int mx = w*GetMouseX()/GetScreenWidth(),my = h*GetMouseY()/GetScreenHeight();
		float s = GetMouseWheelMove();
		cursor.vx = mx - cursor.x;
		cursor.vy = my - cursor.y;
		cursor.s = s;
		cursor.x = mx;
		cursor.y = my;
		if(held_idx != -1){
			if(IsMouseButtonPressed(1)){
				n_panes += 1;
				p_ptr = realloc(p_ptr, n_panes*sizeof(Pane));
				*(p_ptr + (n_panes - 1)) = *init_sub_pane(p_ptr+held_idx, 0, 0, 0.5, 0.5,c_fun_from_pix,n_panes-1,sub_pix,1,1);

			}

			Rectangle ol_r = (p_ptr+held_idx)->r;		
			move_pane(p_ptr+held_idx, &cursor);
	//		move_childs(p_pPane *ptr, &cursor,n_panes,held_idx);;
					
			upd_childs(p_ptr, n_panes, held_idx,ol_r.width,ol_r.height,ol_r.x,ol_r.y);}

	}
	UnloadTexture(tex);
	CloseWindow();
	return 0;
}

