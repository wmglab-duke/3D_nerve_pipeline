#!/usr/bin/env python3.7

"""
The copyrights of this software are owned by Duke University.
Please refer to the LICENSE.txt and README.txt files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent
"""

# builtin
from typing import List, Tuple

# packages
import pymunk.pygame_util
import numpy as np
from shapely.geometry import LineString, Point

# ascent
from src.core import Trace, Slide
from src.utils import Exceptionable, SetupMode, ReshapeNerveMode, Config


class Deformable(Exceptionable):

    def __init__(self, exception_config: list, boundary_start: Trace, boundary_end: Trace, contents: List[Trace]):
        """
        :param exception_config: pre-loaded data
        :param boundary_start: original start trace
        :param boundary_end: end trace
        :param contents: list of traces assumed to be within boundary start, not required to be within boundary end.
        Assumed boundary end will be able to hold all contents.
        
        """

        # init superclass
        Exceptionable.__init__(self, SetupMode.OLD, exception_config)

        # initialize instance variables
        self.start = boundary_start
        self.end = boundary_end
        self.contents = contents

    def deform(self,
               morph_count: int = 100,
               morph_index_step: int = 10,
               render: bool = True,
               minimum_distance: float = 0.0,
               ratio: float = None, comp = False) -> Tuple[List[tuple], List[float]]:
        """
        :param ratio:
        :param morph_count: number of incremental traces including the start and end of boundary
        :param morph_index_step: steps between loops of updating outer boundary, i.e. 1 is every loop,
        2 is every other loop...
        :param render: True if you care to see it happen... makes this method WAY slower
        :param minimum_distance: separation between original inputs
        :return: tuple of a list of total movement vectors and total angle rotated for each fascicle
        """

        # copy the "contents" so multiple deformations are possible
        contents = [trace.deepcopy() for trace in self.contents]

        bounds = self.start.polygon().bounds
        width = int(1.5 * (bounds[2] - bounds[0]))
        height = int(1.5 * (bounds[3] - bounds[1]))

        # offset all the traces to provide for an effective minimum distance for original fascicles
        for trace in contents:
            trace.offset(distance=minimum_distance / 2.0)

        # initialize drawing vars, regardless of whether or not actually rendering
        # these have been moved below (if render...)
        drawing_screen = options = display_dimensions = screen = None

        # initialize the physics space (gravity is 0)
        space = pymunk.Space()

        # referencing the deform_steps method below
        if not comp: 
            morph_steps = [step.pymunk_segments(space) for step in Deformable.deform_steps(self.start,
                                                                                       self.end,
                                                                                       morph_count,
                                                                                       ratio)]
        else: 
            morph_steps = [step.pymunk_segments(space) for step in Deformable.deform_steps_complex(self.start,
                                                                                       self.end,
                                                                                       morph_count,
                                                                                       ratio)]
        # draw the deformation
        if render:
            # packages
            import pygame
            from pygame.locals import KEYDOWN, K_ESCAPE, QUIT, K_SPACE
            from pygame.colordict import THECOLORS

            width = int(1.5 * (bounds[2] - bounds[0]))
            height = int(1.5 * (bounds[3] - bounds[1]))

            aspect = 2

            display_dimensions = int(800 * aspect), 800
            drawing_screen = pygame.display.set_mode(display_dimensions)

            # pygame.init()
            # drawing_screen = pygame.Surface((width, height))
            options = pymunk.pygame_util.DrawOptions(drawing_screen)
            options.shape_outline_color = (0, 0, 0, 255)
            options.shape_static_color = (0, 0, 0, 255)

        # init vector of start positions
        start_positions: List[np.ndarray] = []
        start_rotations: List[float] = []

        # add fascicles bodies to space
        fascicles = [trace.pymunk_poly() for trace in contents]
        for body, shape in fascicles:
            shape.elasticity = 0.0
            space.add(body, shape)
            start_positions.append(np.array(body.position))
            start_rotations.append(body.angle)

        def add_boundary():
            for seg in morph_step:
                seg.elasticity = 0.0
                seg.group = 1
            space.add(*morph_step)

        def step_physics(space: pymunk.Space, count: int):
            dt = 1.0 / 60.0
            for x in range(count):
                space.step(dt)

        running = True
        loop_count = morph_index_step
        morph_index = 0
        morph_step = morph_steps[morph_index]
        add_boundary()

        while running:
            # if the loop count is divisible by the index step, update morph
            if render:
                for event in pygame.event.get():
                    if event.type == QUIT:
                        running = False
                    elif event.type == KEYDOWN and event.key == K_ESCAPE:
                        running = False
                    elif event.type == KEYDOWN and event.key == K_SPACE:
                        pass

            if loop_count % morph_index_step == 0:
                # print('PRINT PRINT PRINT')
                space.remove(*morph_step)
                morph_index += 1
                Deformable.printProgressBar(morph_index,
                                            len(morph_steps),
                                            prefix='\t\tdeforming',
                                            suffix='complete',
                                            length=50)
                # print('\tmorph step {} of {}'.format(morph_index, len(morph_steps)))

                if morph_index == len(morph_steps):
                    running = False
                else:
                    morph_step = morph_steps[morph_index]
                    add_boundary()

            # update physics
            step_physics(space, 3)
            
            # draw screen
            if render:
                drawing_screen.fill(THECOLORS["white"])

                # draw the actual screen
                space.debug_draw(options)

                temp_surf = drawing_screen.copy()
                # drawing_screen.fill((0,0,0))  # here, you can fill the screen with whatever you want to take the place of what was there before
                # drawing_screen.blit(temp_surf, (width/2,-height/2))

                # pygame.display.update()

                # drawing_screen = pygame.transform.scale(drawing_screen, display_dimensions)#, screen)
                pygame.display.flip()
                pygame.display.set_caption('nerve morph step {} of {}'.format(morph_index, len(morph_steps)))

            loop_count += 1

        step_physics(space, 500)

        # get end positions
        end_positions: List[np.ndarray] = []
        end_rotations: List[float] = []
        for body, _ in fascicles:
            end_positions.append(np.array(body.position))
            end_rotations.append(body.angle)

        # return total movement vectors (dx, dy)
        movements = [tuple(end - start) for start, end in zip(start_positions, end_positions)]
        rotations = [end - start for start, end in zip(start_rotations, end_rotations)]
        return movements, rotations
    
    def fiber_deform(self,points,start,
               morph_count: int = 100,
               morph_index_step: int = 10,
               render: bool = True,
               minimum_distance: float = 0.0,
               ratio: float = None, comp = False) -> Tuple[List[tuple], List[float]]:
        """
        :param ratio:
        :param morph_count: number of incremental traces including the start and end of boundary
        :param morph_index_step: steps between loops of updating outer boundary, i.e. 1 is every loop,
        2 is every other loop...
        :param render: True if you care to see it happen... makes this method WAY slower
        :param minimum_distance: separation between original inputs
        :return: tuple of a list of total movement vectors and total angle rotated for each fascicle
        """
        from pymunk import Body, Circle
        from pymunk.query_info import PointQueryInfo
        from pymunk.vec2d import Vec2d
        import matplotlib.pyplot as plt
        plt.figure()
        def plotty(trace,bodies):
            trace.plot()
            for point in bodies:
                plt.scatter(point.position[0],point.position[1])
        def pymunk_body(p):
            bd = Body(mass=1,moment=1)
            bd.position = Vec2d(p[0],p[1])
            return bd
        
        def add_forces(bodies, coeff):
            def force_calc(body,bodies):
                final_vector = Vec2d(0,0)
                for b in bodies:
                    if b!=body:
                        thisvec = body.position-b.position
                        thisvec*=1/thisvec.get_length_sqrd()
                        final_vector+=thisvec
                return final_vector
            for i,body in enumerate(bodies):
                fvec = force_calc(body,bodies)*coeff
                body.apply_force_at_local_point(tuple([fvec[0],fvec[1]]),tuple([0,0]))
        def add_boundary_force(bodies,boundary,boundcoeff):
            for body in bodies:
                __,pq = boundary.point_query(body.position)
                fvec = body.position-pq.point
                fvec*=1/fvec.get_length_sqrd()
                fvec*= boundcoeff
                body.apply_force_at_local_point(tuple([fvec[0],fvec[1]]),tuple([0,0]))
        def add_all_boundary_force(bodies,trace,boundcoeff):
            def force_calc(body,trace):
                final_vector = Vec2d(0,0)
                for p in trace.points:
                    thisvec = body.position-Vec2d(p[0],p[1])
                    thisvec= thisvec*float(1/thisvec.get_length_sqrd())
                    final_vector+=thisvec
                return final_vector
            for i,body in enumerate(bodies):
                fvec = force_calc(body,trace)*coeff
                body.apply_force_at_local_point(tuple([fvec[0],fvec[1]]),tuple([0,0]))
                
        bounds = self.start.polygon().bounds
        width = int(1.5 * (bounds[2] - bounds[0]))
        height = int(1.5 * (bounds[3] - bounds[1]))
        
        loop_n = 30

        def add_boundary():
            for seg in morph_step:
                seg.elasticity = 0.0
                seg.group = 1
            space.add(*morph_step)
            
        # initialize drawing vars, regardless of whether or not actually rendering
        # these have been moved below (if render...)
        drawing_screen = options = display_dimensions = screen = None

        # initialize the physics space (gravity is 0)
        space = pymunk.Space()
        space.damping = .9

        # referencing the deform_steps method below
        morph_steps = Deformable.deform_steps_complex(self.start,self.end,morph_count,ratio)
        morph_index = 0
        morph_step = morph_steps[morph_index]

        # init vector of start positions
        start_positions: List[np.ndarray] = []

        # add fascicles bodies to space
        fibers = [pymunk_body(p) for p in points]
        shapes = [Circle(b,radius = 1) for b in fibers]
        for body in fibers:
            space.add(body)
            start_positions.append(np.array(body.position))
        for shape in shapes: space.add(shape)
        
        running = True
        
        coeff = 1
        boundcoeff = 30
        
        while running:
            nbody, boundary = morph_steps[morph_index].pymunk_poly()
            space.add(boundary)
            for i in range(loop_n):
                add_forces(fibers,coeff)
                # add_boundary_force(fibers, boundary, boundcoeff)
                add_all_boundary_force(fibers, start, boundcoeff)
                # update physics
                space.step(1)
                plt.figure()
                plotty(morph_steps[morph_index],fibers)
            space.remove(*morph_step)
            morph_index += 1
            if morph_index == len(morph_steps):
                running = False
            else:
                morph_step = morph_steps[morph_index]
                add_boundary()
            #flag, dont forget to check for multipoint intersections on the ray that detects the boundaries
            #flag, need to find nearest point on boundary and apply a repulsive force from it
        # get end positions
        end_positions: List[np.ndarray] = []
        for body in fibers:
            end_positions.append(np.array(body.position))

        # return total movement vectors (dx, dy)
        movements = [tuple(end - start) for start, end in zip(start_positions, end_positions)]
        return movements

    @staticmethod
    def deform_steps(start: Trace, end: Trace, count: int = 2, deform_ratio: float = 1.0, slide: Slide = None) \
            -> List[Trace]:

        # Find point along old_nerve that is closest to major axis of best fit ellipse
        (x, y), (a, b), angle = start.ellipse()  # returns degrees

        angle *= 2 * np.pi / 360  # converts to radians

        ray = LineString([(x, y), (x + (2 * a * np.cos(angle)), y + (2 * a * np.sin(angle)))])

        start_intersection = ray.intersection(start.polygon().boundary)
        end_intersection = ray.intersection(end.polygon().boundary)

        start_distances = [Point(point[:2]).distance(start_intersection) for point in start.points]
        end_distances = [Point(point[:2]).distance(end_intersection) for point in end.points]

        # Use point on major axis as the first point on old_nerve, find closest point on new_nerve
        # and assign as first point

        # Sweep CW assigning consecutive points
        start_initial_index = np.where(np.array(start_distances == np.min(start_distances)))[0][0]
        end_initial_index = np.where(np.array(end_distances == np.min(end_distances)))[0][0]

        start_i = start_initial_index
        end_i = end_initial_index
        associated_points = [(0, 0)] * len(start.points)

        # Match pre and post-deformation nerve points. Increment one, decrement the other
        # to match rotation of trace points.
        while True:
            associated_points[start_i] = (start.points[start_i], end.points[end_i])

            # map slide.orientation_point_index
            if (slide is not None) and \
                    (slide.orientation_point_index is not None) and \
                    (slide.orientation_point_index == start_i):
                slide.orientation_point_index = end_i

            if start_i == len(start.points) - 1:
                start_i = 0
            else:
                start_i += 1

            if end_i == 0:
                end_i = len(end.points) - 1
            else:
                end_i -= 1

            if start_i == start_initial_index and end_i == end_initial_index:
                break

        # Find vector between old_nerve and new_nerve associated points
        vectors = [second - first for first, second in associated_points]

        # Save incremental steps of nerve
        ratios = np.linspace(0, 1, count)
        traces = []
        for ratio in ratios:
            trace = start.deepcopy()
            for i, point in enumerate(trace.points):
                point += vectors[i] * ratio
            traces.append(trace)

        return traces[:int((deform_ratio if deform_ratio is not None else 1) * count)]

    @staticmethod
    def deform_steps_complex(start: Trace, end: Trace, count: int = 2, deform_ratio: float = 1.0, slide: Slide = None, plot = False) \
            -> List[Trace]:
        from shapely.geometry import LinearRing
        from shapely.ops import split
        import sys
        import matplotlib.pyplot as plt
        def plotty(in1,in2):
            x1 = list(zip(*[x for x in in1.coords]))[0]
            y1 = list(zip(*[x for x in in1.coords]))[1]
            x2 = list(zip(*[x for x in in2.coords]))[0]
            y2 = list(zip(*[x for x in in2.coords]))[1]
            plt.plot(x1,y1)
            plt.plot(x2,y2)
        start = start.deepcopy()
        end = end.deepcopy()
        shift = list(np.array([0,0]) - np.array(start.centroid())) + [0]
        start.shift(shift)
        shift = list(np.array([0,0]) - np.array(end.centroid())) + [0]
        end.shift(shift)
        startline = LineString([tuple(point) for point in start.points[:, :2]]+[tuple(start.points[0, :2])])
        endline = LineString([tuple(point) for point in end.points[:, :2]]+[tuple(end.points[0, :2])])
        # plotty(startline,endline)
        # plt.scatter(endline.coords[0][0],endline.coords[0][1])

        if not(startline.is_ring and endline.is_ring): sys.exit('whoopsie')
        # Find point along old_nerve that is closest to major axis of best fit ellipse
        (x, y), (a, b), angle = start.ellipse()  # returns degrees

        angle *= 2 * np.pi / 360  # converts to radians

        ray = LineString([(x, y), (x + (2 * a * np.cos(angle)), y + (2 * a * np.sin(angle)))])
        
        start_intersection = ray.intersection(start.polygon().boundary)
        #assert that is single point because currently there is no logic to handle concave shapes
        #however all getting this working would involve is making the splitting logic below only split at the point that 
        #ray intersects the boundary which is closest to the centroid
        assert isinstance(start_intersection,Point)
        end_intersection = ray.intersection(end.polygon().boundary)
        #assert again (I read that in the the voice of the cha cha slide guy)
        assert isinstance(end_intersection,Point)
        #need to add check to make sure intersection is only one point
        startsplit = split(startline,ray)
        endsplit = split(endline,ray)
        # plotty(startsplit[1],ray)
        # plt.scatter(startsplit[0].coords[0][0],startsplit[0].coords[0][1])
        startfix = LineString(list(startsplit[1].coords)+list(startsplit[0].coords))
        endfix = LineString(list(endsplit[1].coords)+list(endsplit[0].coords))
        # plotty(startline,ray)
        # plt.scatter(endfix.coords[0][0],endfix.coords[0][1])
        # plt.scatter(endfix.coords[-1][0],endfix.coords[-1][1])

        interpcount = 500
        # Find vector between old_nerve and new_nerve associated points
        ns = [np.array(startfix.interpolate(d,normalized = True).coords[0]) for d in np.linspace(0,1,num=interpcount)]
        ns = [np.append(x,0) for x in ns]
        ne = [np.array(endfix.interpolate(d,normalized = True).coords[0]) for d in np.linspace(0,1,num=interpcount)]
        ne = [np.append(x,0) for x in ne] 
        vectors = [ne[i]-ns[i] for i in range(len(ne))]
            
        # Save incremental steps of nerve
        ratios = np.linspace(0, 1, count)
        traces = []
        for ratio in ratios:
            trace = Trace(ns,Config.EXCEPTIONS.value)
            for i, point in enumerate(trace.points):
                point += vectors[i] * ratio
            traces.append(trace)
        # traces[0].plot()
   
        # start.plot()
        # traces[-1].plot()
        # end.plot()
        if plot:
            for trace in traces:
                plt.figure()
                trace.plot()
        return traces[:int((deform_ratio if deform_ratio is not None else 1) * count)]

    @staticmethod
    def from_slide(slide: Slide, mode: ReshapeNerveMode, sep_nerve: float = None) -> 'Deformable':
        # method in slide will pull out each trace and add to a list of contents, go through traces and build polygons

        # exception configuration data
        exception_config_data = slide.configs[Config.EXCEPTIONS.value]

        bounds = slide.nerve.polygon().bounds
        width = int(1.5 * (bounds[2] - bounds[0])) / 2
        height = int(1.5 * (bounds[3] - bounds[1])) / 2

        slide.move_center(np.array([1.5 * width, 1.5 * height]))

        # settle the inners
        # for fascicle in slide.fascicles:
        #     inners = [inner.deepcopy() for inner in fascicle]
        #     def_tmp = Deformable(exception_config_data, fascicle.outer, fascicle.outer, inners)
        #     movements, rotations = def_tmp.deform(morph_count=36,
        #                                           render=False,
        #                                           minimum_distance=10)
        #     for move, angle, inner in zip(movements, rotations, fascicle.inners):
        #         inner.shift(list(move) + [0])
        #         inner.rotate(angle)

        # get start boundary
        boundary_start = slide.nerve.deepcopy()

        # get end boundary
        if sep_nerve is None:
            sep_nerve = 0.0
        boundary_end = slide.reshaped_nerve(mode, buffer=sep_nerve).deepcopy()

        # get contents
        contents = [fascicle.outer.deepcopy() for fascicle in slide.fascicles]

        slide.move_center(np.array([0, 0]))

        # return new object
        return Deformable(exception_config_data, boundary_start, boundary_end, contents)

    # copied from https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    @staticmethod
    def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█'):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
        """
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print('\r{} |{}| {}% {}'.format(prefix, bar, percent, suffix), end='')
        # Print New Line on Complete
        if iteration == total:
            print()
