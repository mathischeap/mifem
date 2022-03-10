import random

num_containers = 10
num_objects = 500

containers = [[] for _ in range(num_containers)]
objects = (1 for _ in range(num_objects))
# we don't need a list of these, just have to iterate over it, so this is a genexp

for object in objects:
    random.choice(containers).append(object)


LIST = list()
for c in containers:
    LIST.append(len(c))

print(LIST)
print(sum(LIST))


a = random.randint(0,2)

print(a)